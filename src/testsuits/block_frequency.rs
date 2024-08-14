use super::util::*;
use crate::{Sample, TestResult};

/// 块内频数检测
pub(crate) fn block_frequency(sample: &Sample, m: i32) -> TestResult {
    let mut v = 0.0;
    let m = m as usize;

    let mut i = 0;
    while i <= sample.e.len() - m {
        let pi = popcount(&sample.e[i..i + m]) as f64 / m as f64 - 0.5;
        v += pi * pi;
        i += m;
    }

    let n = (sample.e.len() / m) as f64;
    let pv = igamc(n / 2.0, v * 2.0 * m as f64);
    TestResult {
        pv,
        qv: pv,
    }
}



