use super::util::*;
use crate::{Sample, TestResult};

/// 游程总数检测
pub(crate) fn runs(sample: &Sample) -> TestResult {
    let e = &sample.e;
    let n = e.len();
    let mut v = 1u64;
    for i in 0..n - 1 {
        v += (e[i] ^ e[i + 1]) as u64;
    }
    let pi = sample.pop as f64 / n as f64;
    let t = 2.0 * pi * (1.0 - pi);
    let pv = erfc(abs(v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)));
    let qv = erfc((v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)))/ 2.0;
    TestResult {
        pv1: pv,
        qv1: qv,
        pv2: None,
        qv2: None,
    }
}