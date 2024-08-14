use super::util::*;
use crate::{Sample, TestResult};

/// 游程分布检测
pub(crate) fn runs_distribution(sample: &Sample) -> TestResult {
    // n          k
    // 100      : 2
    // 1000     : 5
    // 10000    : 8
    // 100000   : 12
    // 1000000  : 15
    let n = sample.e.len();
    let mut ei;
    let mut k = 0;
    for i in 1..(n + 1) {
        ei = (n - i + 3) as f64 / powi(2.0, i as i32 + 2);
        if ei >= 5.0 {
            k = if i > k { i } else { k };
        }
        if ei < 1.0 {
            break;
        }
    }

    // bi[0] and gi[0] are dummy.
    let mut bi = vec![0; k+1];
    let mut gi = vec![0; k+1];

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
    for i in 1..k {
        ei = t / powi(2.0, i as i32 + 1);
        v += (powi(bi[i] as f64 - ei, 2) + powi(gi[i] as f64 - ei, 2)) / ei;
    }
    ei = t / powi(2.0, k as i32);
    v += (powi(bi[k] as f64 - ei, 2) + powi(gi[k] as f64 - ei, 2)) / ei;

    let pv = igamc((k - 1) as f64, v / 2.0);
    TestResult {
        pv,
        qv: pv,
    }

}