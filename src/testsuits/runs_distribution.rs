use super::util::*;
use crate::{Sample, TestResult};

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
    let mut bi = vec![0; k + 1];
    let mut gi = vec![0; k + 1];

    // current_run > 0 means the 1's run in the current
    // current_run < 0 means the 0's run in the current
    let mut current_run = 2 * sample.get_bit_unchecked(0) as i32 - 1;
    for i in 1..n {
        if sample.get_bit_unchecked(i) == 1 {
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

#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

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
        for nbits in 9..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(
                runs_distribution_u64(&sample),
                runs_distribution_epsilon(&sample)
            );
        }
    }

    // #[test]
    // fn test_ei() {
    //     let n = 100000;
    //     println!("{}:{}", n, ei(n));
    // }
}
