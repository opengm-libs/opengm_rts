use super::util::*;
use crate::{Sample, TestResult};

/// 自相关检测
pub(crate) fn autocorrelation(sample: &Sample, d: i32) -> TestResult{
    let epsilon = &sample.e;
    let n = epsilon.len();
    let d = d as usize;

    let mut sum = 0i64;
	for i in 0..(n-d){
		sum += (epsilon[i] ^ epsilon[i+d]) as i64;
	}

    let v = (2.0 * sum as f64 - (n-d) as f64) / sqrt((n-d) as f64);
	let pv = erfc(abs(v)/SQRT2);
	let qv = erfc(v/SQRT2)/2.0;

    TestResult {
        pv1: pv,
        qv1: qv,
        pv2: None,
        qv2: None,
    }
}