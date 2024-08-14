use super::util::*;
use crate::{Sample, TestResult};

/// 累加和检测
pub(crate) fn cumulative_sums_forward(sample: &Sample) -> TestResult {
    let pv = cumulative_sums_inner(&sample.e, true);
    TestResult {
        pv,
        qv: pv,
    }
}

pub(crate) fn cumulative_sums_backward(sample: &Sample) -> TestResult {
    let pv = cumulative_sums_inner(&sample.e, false);
    TestResult {
        pv,
        qv: pv,
    }
}

fn cumulative_sums_inner(e: &[u8], forward: bool) -> f64 {
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
