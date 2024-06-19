use super::util::*;
use crate::Sample;

fn psi2(e: &[u8], m: i32) -> f64 {
    let n = e.len();
    // m is at most 7 in GM/T 0005-2021 when n = 100 000 000.
    const MAX_M: usize = 1 << 10;
    if m <= 0 {
        return 0.0;
    }

    let length = 1 << m;
    let mask = length - 1;
    let mut bucket = [0u64; MAX_M];

    let mut k = 0;
    let mut i: usize = 0;
    while i < m as usize {
        k = (k << 1) | (e[i] as usize);
        i += 1;
    }
    bucket[k] += 1;

    for _ in 1..n {
        k = (k << 1) | (e[i] as usize);
        k &= mask;
        bucket[k] += 1;
        i += 1;
        if i == n {
            i = 0;
        }
    }

    // Note that a uint sum may overflow.
    let mut sum = 0.0;
    for i in 0..length {
        sum += (bucket[i] * bucket[i]) as f64;
    }
    sum / (n as f64) * (1 << m) as f64 - n as f64
}

// 4. 重叠子序列检测
// 2 <= m <= 10
pub(crate) fn serial(sample: &Sample, m: i32) -> (f64, f64) {
    if m < 2 || m > 10 {
        panic!("serial: m invalid\n");
    }

    let p0 = psi2(&sample.e, m);
    let p1 = psi2(&sample.e, m - 1);
    let p2 = psi2(&sample.e, m - 2);
    let del1 = p0 - p1;
    let del2 = p0 - 2.0 * p1 + p2;

    let pvalue1 = igamc(powi(2.0, m - 2), del1 / 2.0);
    let pvalue2 = igamc(powi(2.0, m - 3), del2 / 2.0);

    (pvalue1, pvalue2)
}