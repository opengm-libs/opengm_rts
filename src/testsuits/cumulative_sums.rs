use super::util::*;
use crate::{get_bit_unchecked, Sample, TestResult};


pub(crate) fn cumulative_sums_forward(sample: &Sample) -> TestResult {
    cumulative_sums_forward_u64(sample)
}

pub(crate) fn cumulative_sums_backward(sample: &Sample) -> TestResult {
    cumulative_sums_backward_u64(sample)
}

////////////////////////////////////////////////////////////////

/// 累加和检测
#[cfg(test)]
pub(crate) fn cumulative_sums_forward_epsilon(sample: &Sample) -> TestResult {
    let pv = cumulative_sums_inner_epsilon(&sample.e, true);
    TestResult {
        pv,
        qv: pv,
    }
}

#[cfg(test)]
pub(crate) fn cumulative_sums_backward_epsilon(sample: &Sample) -> TestResult {
    let pv = cumulative_sums_inner_epsilon(&sample.e, false);
    TestResult {
        pv,
        qv: pv,
    }
}

fn cumulative_sums_inner_epsilon(e: &[u8], forward: bool) -> f64 {
    let n = e.len() as i32;
    let mut S = 0;
    let mut sup = 0;
    let mut inf = 0;
    let mut z = 0;
  
    if forward{   
        for k in e {
            S += 2 * (*k as i32) - 1;
            if S > sup {
                sup += 1;
            }
            if S < inf {
                inf -= 1;
            }
            z = if sup > -inf { sup } else { -inf };
        }
    }else{
        for k in e.iter().rev() {
            S += 2 * (*k as i32) - 1;
            if S > sup {
                sup += 1;
            }
            if S < inf {
                inf -= 1;
            }
            z = if sup > -inf { sup } else { -inf };
        }
    }

    let mut sum1 = 0.0;
    let mut k = (-n / z + 1) / 4;
    while k <= (n / z - 1) / 4 {
        sum1 += normal(((4 * k + 1) * z) as f64 / sqrt(n as f64));
        sum1 -= normal(((4 * k - 1) * z) as f64 / sqrt(n as f64));
        k += 1;
    }

    let mut sum2 = 0.0;
    let mut k = (-n / z - 3) / 4;
    while k <= (n / z - 1) / 4 {
        sum2 += normal(((4 * k + 3) * z) as f64 / sqrt(n as f64));
        sum2 -= normal(((4 * k + 1) * z) as f64 / sqrt(n as f64));
        k += 1;
    }

    return 1.0 - sum1 + sum2;
}

////////////////////////////////////////////////////////////////

/// 累加和检测
pub(crate) fn cumulative_sums_forward_u64(sample: &Sample) -> TestResult {
    let pv = cumulative_sums_inner_u64(&sample.b64, sample.len(), true);
    TestResult {
        pv,
        qv: pv,
    }
}

pub(crate) fn cumulative_sums_backward_u64(sample: &Sample) -> TestResult {
    let pv = cumulative_sums_inner_u64(&sample.b64, sample.len(), false);
    TestResult {
        pv,
        qv: pv,
    }
}
fn cumulative_sums_inner_u64(b64: &[u64], n:usize, forward: bool) -> f64 {
    let n = n as i32;
    let mut S = 0;
    let mut sup = 0;
    let mut inf = 0;
    let mut z = 0;
  
    if forward{   
        for i in 0..n {
            let k = get_bit_unchecked(b64, i as usize) as i32;
            S += 2 * k - 1;
            if S > sup {
                sup += 1;
            }
            if S < inf {
                inf -= 1;
            }
            z = if sup > -inf { sup } else { -inf };
        }
    }else{
        for i in (0..n).rev() {
            let k = get_bit_unchecked(b64, i as usize) as i32;
            S += 2 * k - 1;
            if S > sup {
                sup += 1;
            }
            if S < inf {
                inf -= 1;
            }
            z = if sup > -inf { sup } else { -inf };
        }
    }

    let mut sum1 = 0.0;
    let mut k = (-n / z + 1) / 4;
    while k <= (n / z - 1) / 4 {
        sum1 += normal(((4 * k + 1) * z) as f64 / sqrt(n as f64));
        sum1 -= normal(((4 * k - 1) * z) as f64 / sqrt(n as f64));
        k += 1;
    }

    let mut sum2 = 0.0;
    let mut k = (-n / z - 3) / 4;
    while k <= (n / z - 1) / 4 {
        sum2 += normal(((4 * k + 3) * z) as f64 / sqrt(n as f64));
        sum2 -= normal(((4 * k + 1) * z) as f64 / sqrt(n as f64));
        k += 1;
    }

    return 1.0 - sum1 + sum2;
}



#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::CumulativeSumsForward);
        let sample: Sample = tv.0.into();
        assert_eq!(cumulative_sums_forward_u64(&sample), tv.2);
        assert_eq!(cumulative_sums_forward_epsilon(&sample), tv.2);

        let tv = get_test_vec(crate::TestFuncs::CumulativeSumsBackward);
        let sample: Sample = tv.0.into();
        assert_eq!(cumulative_sums_backward_epsilon(&sample), tv.2);
        assert_eq!(cumulative_sums_backward_u64(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::CumulativeSumsForward);
        let sample: Sample = E.into();
        assert_eq!(tv.1, cumulative_sums_forward_epsilon(&sample));
        assert_eq!(tv.1, cumulative_sums_forward_u64(&sample));

        let tv = get_test_vec_e(crate::TestFuncs::CumulativeSumsBackward);
        let sample: Sample = E.into();
        assert_eq!(tv.1, cumulative_sums_backward_epsilon(&sample));
        assert_eq!(tv.1, cumulative_sums_backward_u64(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 128/8..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(cumulative_sums_backward_epsilon(&sample), cumulative_sums_backward_u64(&sample));
            assert_eq!(cumulative_sums_forward_epsilon(&sample), cumulative_sums_forward_u64(&sample));
        }
    }
}
