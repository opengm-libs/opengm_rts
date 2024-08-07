use super::util::*;
use crate::{Sample, TestResult};

const L: usize = 7;
const Q: usize = 1280;

/// 通用检测
/// Note that sample must have length greater than 7*1280.
pub(crate) fn universal(sample: &Sample) -> TestResult {
    const E: f64 = 6.1962507;
    const VAR: f64 = 3.125;
    let n = sample.e.len();
    if n < L * Q {
        return TestResult::default();
    }
    let K = n / L - Q;
    let p = 1 << L;
    let mut T = vec![0usize; p];

    let c =
        0.7 - 0.8 / (L as f64) + (4.0 + 32.0 / (L as f64)) * powf(K as f64, -3.0 / L as f64) / 15.0;
    let sigma = c * sqrt(VAR / K as f64);
    let mut sum = 0.0;
    let mut dec_rep;
    for i in 1..=Q {
        dec_rep = 0;
        for j in 0..L {
            dec_rep += sample.e[(i - 1) * L + j] as usize * (1 << (L - 1 - j));
        }
        T[dec_rep] = i;
    }
    for i in (Q + 1)..=(Q + K) {
        dec_rep = 0;
        for j in 0..L {
            dec_rep += sample.e[(i - 1) * L + j] as usize * (1 << (L - 1 - j));
        }
        sum += ln(i as f64 - (T[dec_rep as usize] as f64)) / ln(2.0);
        T[dec_rep] = i;
    }
    let phi = sum / K as f64;
    let V = (phi - E) / sigma;

    let pv = erfc(abs(V) / SQRT2);
    let qv = erfc(V / SQRT2) / 2.0;
    TestResult {
        pv1: pv,
        qv1: qv,
        pv2: None,
        qv2: None,
    }
}
