use super::util::*;
use crate::{get_bit_unchecked, Sample, TestResult};

/// 游程分布检测
pub(crate) fn runs_distribution(sample: &Sample) -> TestResult {
    runs_distribution_u64(sample)
}

#[inline]
fn ei(n: usize) -> usize {
    match n {
        100 => 2,
        1000 => 5,
        10000 => 8,
        100000 => 12,
        1000000 => 15,
        10000000 => 18,
        100000000 => 22,
        _ => {
            let mut e;
            let mut k = 0;
            for i in 1..=n {
                e = (n - i + 3) as f64 / powi(2.0, i as i32 + 2);
                if e >= 5.0 {
                    k = if i > k { i } else { k };
                }
                if e < 1.0 {
                    break;
                }
            }
            k
        }
    }
}

#[cfg(test)]
pub(crate) fn runs_distribution_epsilon(sample: &Sample) -> TestResult {
    let n = sample.e.len();
    let k = ei(n);

    // bi[0] and gi[0] are dummy.
    let mut bi = vec![0; k + 1];
    let mut gi = vec![0; k + 1];

    // current_run > 0 means the 1's run in the current
    // current_run < 0 means the 0's run in the current
    let mut current_run = 2 * sample.e[0] as i32 - 1;
    for i in 1..n {
        if sample.e[i] == 1 {
            if current_run < 0 {
                // in the case 0...01, record the 0's run.
                gi[saturating(-current_run, 0, k as i32) as usize] += 1;
                current_run = 1;
            } else {
                // in the case 0...00
                current_run += 1;
            }
        } else {
            if current_run > 0 {
                bi[saturating(current_run, 0, k as i32) as usize] += 1;
                current_run = -1;
            } else {
                current_run -= 1;
            }
        }
    }
    // Record the last run
    if current_run > 0 {
        bi[saturating(current_run, 0, k as i32) as usize] += 1;
    } else {
        gi[saturating(-current_run, 0, k as i32) as usize] += 1;
    }

    let mut t = 0f64;
    for i in 1..(k + 1) {
        t += (bi[i] + gi[i]) as f64;
    }

    let mut v = 0.0;
    let mut ei;
    for i in 1..k {
        ei = t / powi(2.0, i as i32 + 1);
        v += (powi(bi[i] as f64 - ei, 2) + powi(gi[i] as f64 - ei, 2)) / ei;
    }
    ei = t / powi(2.0, k as i32);
    v += (powi(bi[k] as f64 - ei, 2) + powi(gi[k] as f64 - ei, 2)) / ei;

    let pv = igamc((k - 1) as f64, v / 2.0);
    TestResult { pv, qv: pv }
}

////////////////////////////////////////////////////////////////
///
pub(crate) fn runs_distribution_u64(sample: &Sample) -> TestResult {
    let n = sample.bit_length;
    let k = ei(n);
    // let b64 = sample.b64.as_slice();

    // bi[0] and gi[0] are dummy.
    let mut bg = vec![vec![0; k + 1], vec![0; k + 1]];
    count_runs(&sample.b64, sample.len(), &mut bg, k);
    let bi = &bg[0];
    let gi = &bg[1];

    let mut t = 0f64;
    for i in 1..(k + 1) {
        t += (bi[i] + gi[i]) as f64;
    }

    let mut v = 0.0;
    let mut ei;
    for i in 1..k {
        ei = t / powi(2.0, i as i32 + 1);
        v += (powi(bi[i] as f64 - ei, 2) + powi(gi[i] as f64 - ei, 2)) / ei;
    }
    ei = t / powi(2.0, k as i32);
    v += (powi(bi[k] as f64 - ei, 2) + powi(gi[k] as f64 - ei, 2)) / ei;

    let pv = igamc((k - 1) as f64, v / 2.0);
    TestResult { pv, qv: pv }
}

////////////////////////////////////////////////////////////////
fn count_runs(b64: &[u64], n: usize, runs: &mut Vec<Vec<u64>>, k: usize) {
    let mut tail = 1 ^ ((b64[0] >> 63) as u8);
    let mut tail_run = 0;
    if n % 64 == 0 {
        for b in b64 {
            (tail, tail_run) = count_runs_single(tail, tail_run, *b, runs, k);
        }
    } else {
        for b in b64[..b64.len() - 1].iter() {
            (tail, tail_run) = count_runs_single(tail, tail_run, *b, runs, k);
        }
        let mut x = b64[b64.len() - 1];
        let last_bit = get_bit_unchecked(b64, n - 1);
        let last_run = 64 - n % 64;
        if last_bit == 0 {
            x |= lower_bits_mask_u64(last_run);
        }
        count_runs_single(tail, tail_run, x, runs, k);
        let idx = saturating_ceil(last_run, k);
        runs[1 ^ last_bit as usize][idx] -= 1;
    }
}

fn count_runs_single(
    last_bit: u8,
    last_run: u32,
    x: u64,
    runs: &mut Vec<Vec<u64>>,
    k: usize,
) -> (u8, u32) {
    let res = x as u8 & 1;
    let mut leading;
    let mut x = x;
    let x0 = x & 1;
    if x0 == 0 {
        x = !x;
    }
    let mut n = 0;
    let mut i = last_bit as usize;

    // The first run should be handled separately.
    leading = if last_bit != x0 as u8 {
        x.leading_zeros()
    } else {
        x.leading_ones()
    };

    runs[i][saturating_ceil(leading as usize + last_run as usize, k)] += 1;
    i ^= 1;
    n += leading;
    x = x.wrapping_shl(leading);

    // fix
    runs[last_bit as usize][saturating_ceil(last_run as usize, k)] -= 1;

    while n < 64 {
        // Note that one of leading_xxxs is 0.
        leading = x.leading_ones() + x.leading_zeros();
        let idx = saturating_ceil(leading as usize, k);
        runs[i][idx] += 1;
        i ^= 1;
        n += leading;
        // be careful: if x = 0 or 0xff..fff, then << 64 will panic.
        x = x.wrapping_shl(leading);
    }
    (res, leading)
}

#[cfg(test)]
mod tests {
    use super::{super::tests::*, *};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_runs_u64() {
        let k = 10;
        let mut bg = vec![vec![1; k + 1], vec![1; k + 1]];
        let bit = count_runs_single(0, 0, 0, &mut bg, 10);
        println!("{}", bit.0);
    }

    #[test]
    fn test_runs() {
        let k = 10;
        let sample: Sample = E.into();
        let mut bg = vec![vec![0; k + 1], vec![0; k + 1]];
        count_runs(&sample.b64, sample.len(), &mut bg, k);
    }

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::RunsDistribution);
        let sample: Sample = tv.0.into();
        assert_eq!(runs_distribution_epsilon(&sample), tv.2);
        assert_eq!(runs_distribution_u64(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv1 = get_test_vec_e(crate::TestFuncs::RunsDistribution);
        let sample: Sample = E.into();
        assert_eq!(tv1.1, runs_distribution_epsilon(&sample));
        assert_eq!(tv1.1, runs_distribution_u64(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 9..2000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(
                runs_distribution_u64(&sample),
                runs_distribution_epsilon(&sample)
            );
        }
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::{test_data::E, Sample};
    use test::Bencher;
    #[bench]
    fn bench_runs(b: &mut Bencher) {
        let k = 10;
        let sample: Sample = E.into();
        let mut bg = vec![vec![0; k + 1], vec![0; k + 1]];
        // m= 952,479.20
        b.iter(|| {
            test::black_box(count_runs(&sample.b64, sample.len(), &mut bg, k));
        });
    }
}
