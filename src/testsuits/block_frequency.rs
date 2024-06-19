use super::util::*;
use crate::Sample;

// 2. 块内频数检测
pub(crate) fn block_frequency(sample: &Sample, m: i32) -> f64 {
    let mut v = 0.0;
    let m = m as usize;

    let mut i = 0;
    while i <= sample.e.len() - m {
        let pi = popcount(&sample.e[i..i + m]) as f64 / m as f64 - 0.5;
        v += pi * pi;
        i += m;
    }

    let n = (sample.e.len() / m) as f64;
    igamc(n / 2.0, v * 2.0 * m as f64)
}



