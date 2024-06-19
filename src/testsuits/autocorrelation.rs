use super::util::*;
use crate::Sample;

pub(crate) fn autocorrelation(sample: &Sample, d: i32) -> f64{
    let epsilon = &sample.e;
    let n = epsilon.len();
    let d = d as usize;

    let mut sum = 0i64;
	for i in 0..(n-d){
		sum += (epsilon[i] ^ epsilon[i+d]) as i64;
	}

    let v = (2.0 * sum as f64 - (n-d) as f64) / sqrt((n-d) as f64);
	return  erfc(abs(v)/SQRT2);
}