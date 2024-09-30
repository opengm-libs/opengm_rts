use super::util::*;
use crate::{Sample, TestResult};

/// 游程总数检测
pub(crate) fn runs(sample: &Sample) -> TestResult {
    runs_u64(sample)
}

#[cfg(test)]
pub(crate) fn runs_epsilon(sample: &Sample) -> TestResult {
    let e = &sample.e;
    let n = e.len();
    let mut v = 1u64;
    for i in 0..n - 1 {
        v += (e[i] ^ e[i + 1]) as u64;
    }
    let pi = sample.pop as f64 / n as f64;
    let t = 2.0 * pi * (1.0 - pi);
    let pv = erfc(abs(v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)));
    let qv = erfc((v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64))) / 2.0;
    TestResult { pv, qv }
}

pub(crate) fn runs_u64(sample: &Sample) -> TestResult {
    let b64 = &sample.b64;
    let n = sample.bit_length;
    let mut v = 1u64;
    for i in 0..b64.len()-1{
        let a = b64[i] ^((b64[i] << 1) | (b64[i+1] >> 63));
        v += a.count_ones() as u64;
    }
    let mut last = b64[b64.len()-1];
    last ^= last << 1;
    let tail_length = n%64;
    // clear the last bit
    last &= !(1<<((64-tail_length)%64));
    v += last.count_ones() as u64;

    let pi = sample.pop as f64 / n as f64;
    let t = 2.0 * pi * (1.0 - pi);
    let pv = erfc(abs(v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)));
    let qv = erfc((v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64))) / 2.0;
    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Runs);
        let sample: Sample = tv.0.into();
        assert_eq!(runs_u64(&sample), tv.2);
        assert_eq!(runs_epsilon(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv1 = get_test_vec_e(crate::TestFuncs::Runs);
        let sample: Sample = E.into();
        assert_eq!(tv1.1, runs_epsilon(&sample));
        assert_eq!(tv1.1, runs_u64(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 9..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(runs_epsilon(&sample), runs_u64(&sample));
        }
    }
}


#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::{test_data::E, Sample};
    use test::Bencher;

    #[bench]
    fn bench_test(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 10,432.90 ns/iter
        b.iter(|| {
            test::black_box(runs_u64(&sample));
        });
    }

    #[bench]
    fn bench_test_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 631,445.85 ns/iter
        b.iter(|| {
            test::black_box(runs_epsilon(&sample));
        });
    }
}
