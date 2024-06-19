use super::util::*;
use crate::Sample;

// 1. 频度检测
pub(crate) fn frequency(sample: &Sample) -> f64 {
    let s = sample.pop as i64 * 2 - (sample.e.len() as i64);
    erfc(abs(s as f64) / sqrt((sample.e.len() * 2) as f64))
}