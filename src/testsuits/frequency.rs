use super::util::*;
use crate::{Sample, TestResult};

/// 频度检测
pub(crate) fn frequency(sample: &Sample) -> TestResult {
    let s = sample.pop as i64 * 2 - (sample.e.len() as i64);
    let pv = erfc(abs(s as f64) / sqrt((sample.e.len() * 2) as f64));
    let qv = erfc(s as f64 / sqrt((sample.e.len() * 2) as f64)) / 2.0;
    TestResult {
        pv,
        qv,
    }
}
