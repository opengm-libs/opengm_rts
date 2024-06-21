use super::util::*;
use crate::{Sample, TestResult};

use rustfft::{
    num_complex::{Complex, ComplexFloat},
    FftPlanner,
};

/// 离散傅里叶检测
pub(crate) fn discrete_fourier(sample: &Sample) -> TestResult {
    let e = &sample.e;
    let n = e.len();
    let fft = FftPlanner::<f64>::new().plan_fft_forward(n);
    let mut f = Vec::with_capacity(n);

    for i in 0..n {
        f.push(Complex{ re: 2.0 * e[i] as f64 - 1.0, im: 0.0});
    }
    fft.process(&mut f);
    let t = sqrt(2.995732274 * n as f64);
    let n0 = 0.95 * n as f64 / 2.0;
    let mut n1 = 0;
    for x in &f[0..n / 2 - 1] {
        if x.abs() < t {
            n1 += 1;
        }
    }

    let d = (n1 as f64 - n0) / sqrt(n as f64 * 0.95 * 0.05 / 3.8);
    let pv = erfc(abs(d) / sqrt(2.0));
    let qv = erfc(d/ sqrt(2.0))/2.0;

    TestResult {
        pv1: pv,
        qv1: qv,
        pv2: None,
        qv2: None,
    }
}
