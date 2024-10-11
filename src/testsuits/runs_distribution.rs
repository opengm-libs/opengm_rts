use super::{util::*, super::USE_U8};
use crate::{get_bit_unchecked_u64, Sample, TestResult};

/// 游程分布检测
pub(crate) fn runs_distribution(sample: &Sample) -> TestResult {
    if USE_U8 {
        runs_distribution_u8(sample)
    } else {
        runs_distribution_u64(sample)
    }
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

pub(crate) fn runs_distribution_u64(sample: &Sample) -> TestResult {
    let n = sample.bit_length;
    let k = ei(n);
    // let b64 = sample.b64.as_slice();

    // bi[0] and gi[0] are dummy.
    let mut bg = vec![vec![0; k + 1], vec![0; k + 1]];
    count_runs_u64(&sample.b64, sample.len(), &mut bg, k);
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

fn count_runs_u64(b64: &[u64], n: usize, runs: &mut Vec<Vec<u64>>, k: usize) {
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
        let last_bit = get_bit_unchecked_u64(b64, n - 1);
        let last_run = 64 - n % 64;
        if last_bit == 0 {
            x |= lower_bits_mask_u64(last_run);
        }
        count_runs_single(tail, tail_run, x, runs, k);
        let idx = saturating_ceil(last_run, k);
        runs[1 ^ last_bit as usize][idx] -= 1;
    }
}

////////////////////////////////////////////////////////////////
///

pub(crate) fn runs_distribution_u8(sample: &Sample) -> TestResult {
    let n = sample.bit_length;
    let k = ei(n);

    // bi[0] and gi[0] are dummy.
    let mut bg = vec![vec![0; k + 1], vec![0; k + 1]];
    count_runs_u8(&sample.b, &mut bg, k);
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

fn count_runs_u8(b8: &[u8], runs: &mut Vec<Vec<u64>>, k: usize) {
    // tail init to the flipped bit of epsilon[0].
    let mut last_bit = 1 ^ ((b8[0] >> 7) as u8);
    let mut last_run = 0;

    let full_chunks = b8.len() & (!7);
    for chunk in b8[..full_chunks].chunks_exact(8) {
        let x = u64::from_be_bytes(chunk.try_into().unwrap());
        (last_bit, last_run) =
            count_runs_single(last_bit, last_run, x, runs, k);
    }
    if b8.len() > full_chunks {
        // tail || 00...0
        let mut tail = u64_from_be_slice(&b8[full_chunks..]);
        let tail_bits = (b8.len() - full_chunks) * 8;

        // if tail ends with 0, then x = tail || 11..1
        if (tail >> (64 - tail_bits)) & 1 == 0 {
            tail |= lower_bits_mask_u64(64 - tail_bits);
        }
        count_runs_single(last_bit, last_run, tail, runs, k);

        // fix
        let idx = saturating_ceil(64 - tail_bits, k);
        runs[(tail & 1) as usize][idx] -= 1;
    }
}

////////////////////////////////////////////////////////////////

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
        count_runs_u64(&sample.b64, sample.len(), &mut bg, k);
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
        assert_eq!(tv1.1, runs_distribution_u8(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 13..2000 {
            println!("{:?}", nbits);
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(
                runs_distribution_u64(&sample),
                runs_distribution_epsilon(&sample)
            );
            assert_eq!(
                runs_distribution_u64(&sample),
                runs_distribution_u8(&sample)
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
    fn bench_runs_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m= 4,153,120.85
        b.iter(|| {
            test::black_box(runs_distribution_epsilon(&sample));
        });
    }

    #[bench]
    fn bench_runs_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m= 942,247.90 ns/iter
        b.iter(|| {
            test::black_box(runs_distribution_u64(&sample));
        });
    }

    #[bench]
    fn bench_runs_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m= 976,596.35 ns/iter
        b.iter(|| {
            test::black_box(runs_distribution_u8(&sample));
        });
    }
}
