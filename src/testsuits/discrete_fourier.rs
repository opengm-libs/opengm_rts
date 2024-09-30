use super::util::*;
use crate::{Sample, TestResult};

use rustfft::{
    num_complex::{Complex, ComplexFloat},
    FftPlanner,
};

/// 离散傅里叶检测
pub(crate) fn discrete_fourier(sample: &Sample) -> TestResult {
    discrete_fourier_u64(sample)
}


#[cfg(test)]
pub(crate) fn discrete_fourier_epsilon(sample: &Sample) -> TestResult {
    let e = &sample.e;
    let n = e.len();
    let fft = FftPlanner::<f64>::new().plan_fft_forward(n);
    let mut f = Vec::with_capacity(n);

    for i in 0..n {
        f.push(Complex {
            re: 2.0 * e[i] as f64 - 1.0,
            im: 0.0,
        });
    }
    fft.process(&mut f);

    let t = sqrt(2.995732274 * n as f64);
    let n0 = 0.95 * n as f64 / 2.0;
    let n1 = f[0..=(n / 2 - 1)]
        .iter()
        .map(|x| if x.abs() < t { 1 } else { 0 })
        .sum::<usize>();

    let d = (n1 as f64 - n0) / sqrt(n as f64 * 0.95 * 0.05 / 3.8);
    let pv = erfc(abs(d) / sqrt(2.0));
    let qv = erfc(d / sqrt(2.0)) / 2.0;

    TestResult { pv, qv }
}


#[inline(always)]
fn fft(f: &mut [Complex<f64>]){
    FftPlanner::<f64>::new().plan_fft_forward(f.len()).process(f);
}


/// 离散傅里叶检测
pub(crate) fn discrete_fourier_u64(sample: &Sample) -> TestResult {
    let n = sample.len();
    let mut f = Vec::with_capacity(n);
    for b in sample.b64[..sample.b64.len() - 1].iter() {
        let x = *b;
        for i in (0..64).rev() {
            f.push(Complex {
                re: ((x >> i) & 1) as f64,
                im: 0.0,
            });
        }
    }
    // 0 < tail_length <= 64
    let tail_length = 64 - (64 - sample.len() % 64) % 64;
    let x = sample.b64[sample.b64.len() - 1] >> (64 - tail_length);

    for i in (0..tail_length).rev() {
        f.push(Complex {
            re: ((x >> i) & 1) as f64,
            im: 0.0,
        });
    }

    fft(&mut f);

    let t = sqrt(2.995732274 * n as f64);
    let n0 = 0.95 * n as f64 / 2.0;
    // f0 = sum(x_k) = 2*pop - n
    let mut n1 = if ((2 * sample.pop as i64 - n as i64).abs() as f64) < t {
        1
    } else {
        0
    };
    n1 += f[1..=(n / 2 - 1)]
        .iter()
        .map(|x| if x.abs() < t / 2.0 { 1 } else { 0 })
        .sum::<i64>();

    let d = (n1 as f64 - n0) / sqrt(n as f64 * 0.95 * 0.05 / 3.8);
    let pv = erfc(abs(d) / sqrt(2.0));
    let qv = erfc(d / sqrt(2.0)) / 2.0;

    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use super::{super::tests::*, *};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::DiscreteFourier);
        let sample: Sample = tv.0.into();
        assert_eq!(discrete_fourier_epsilon(&sample), tv.2);
        assert_eq!(discrete_fourier_u64(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::DiscreteFourier);
        let sample: Sample = E.into();
        assert_eq!(tv.1, discrete_fourier_epsilon(&sample));
        assert_eq!(tv.1, discrete_fourier_u64(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 128 / 8..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(
                discrete_fourier_epsilon(&sample),
                discrete_fourier_u64(&sample)
            );
        }
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::{test_data::E, Sample};
    use test::Bencher;

    const K: i32 = 3;

    #[bench]
    fn bench_test(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 30,424,179.20 ns/iter
        b.iter(|| {
            test::black_box(discrete_fourier_epsilon(&sample));
        });
    }

    #[bench]
    fn bench_test_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 28,976,304.10 ns/iter
        b.iter(|| {
            test::black_box(discrete_fourier_u64(&sample));
        });
    }

    #[bench]
    fn bench_test_fft(b: &mut Bencher) {
        let sample: Sample = E.into();
        let e = &sample.e;
        let mut f = Vec::with_capacity(e.len());
        for i in 0..e.len() {
            f.push(Complex {
                re: e[i] as f64,
                im: 0.0,
            });
        }

        // 25,085,366.60 ns/iter
        b.iter(|| {
            let mut ff = f.clone();
            test::black_box(fft(&mut ff));
        });
    }

}
