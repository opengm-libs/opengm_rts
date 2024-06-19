use super::util::*;
use crate::Sample;

const MAX_M :usize = 8;

// 3. 扑克检测
pub(crate) fn poker(sample: &Sample, m: i32) -> f64 {
    let m = m as usize;
    let n = sample.e.len() / m;
    let power: usize = 1 << m;
    // let mut ni = unsafe { Box::<[usize]>::new_zeroed_slice(power).assume_init() };
    let mut ni = vec![0; power];

    for i in 0..n {
        let mut idx = 0_usize;
        for j in 0..m {
            idx = 2 * idx + sample.e[i * m + j] as usize;
        }
        ni[idx] += 1;
    }

    let sum = ni.iter().map(|x| x * x).sum::<usize>();
    let v = (power as f64) / (n as f64) * (sum as f64) - (n as f64);
    let p_value = igamc((power as f64 - 1.0) / 2.0, v / 2.0);
    p_value
}
