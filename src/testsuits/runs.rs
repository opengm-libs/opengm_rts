use super::{util::*, super::USE_U8};
use crate::{Sample, TestResult};

/// 游程总数检测
pub(crate) fn runs(sample: &Sample) -> TestResult {
    if USE_U8 {
        runs_u8(sample)
    } else {
        runs_u64(sample)
    }
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

pub(crate) fn runs_u8(sample: &Sample) -> TestResult {
    if true {
        assert!(sample.len() >= 8);
        let b8 = &sample.b;
        let n = sample.len();

        let mut v = 1u64;

        let mut x = u64::from_be_bytes(b8[..8].try_into().unwrap());
        let mut a;
        let full_chunks = b8.len() & (!7);
        for chunk in b8[8..full_chunks].chunks_exact(8) {
            (a, x) = (x, u64::from_be_bytes(chunk.try_into().unwrap()));
            v += (a ^ ((a << 1) | (x >> 63))).count_ones() as u64;
        }

        if full_chunks == b8.len() {
            let mask = !1;
            x ^= x << 1;
            x &= mask;
            v += x.count_ones() as u64;
        } else {
            // process x
            let y = u64_from_be_slice(&b8[full_chunks..]);
            v += (x ^ ((x << 1) | (y >> 63))).count_ones() as u64;
            let tail_bits = (b8.len() - full_chunks) * 8;
            let mask = (!0) << (64 - tail_bits + 1);
            let y = y ^ (y << 1);
            v += (y & mask).count_ones() as u64;
        }

        let pi = sample.pop as f64 / n as f64;
        let t = 2.0 * pi * (1.0 - pi);
        let pv =
            erfc(abs(v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)));
        let qv =
            erfc((v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64))) / 2.0;
        TestResult { pv, qv }
    } else {
        assert!(sample.len() >= 16);
        let b8 = &sample.b;
        let n = sample.len();

        let mut v = 1u64;

        let mut x = u128::from_be_bytes(b8[..16].try_into().unwrap());
        let mut a;
        let full_chunks = b8.len() & (!15);
        for chunk in b8[16..full_chunks].chunks_exact(16) {
            (a, x) = (x, u128::from_be_bytes(chunk.try_into().unwrap()));
            v += (a ^ ((a << 1) | (x >> 127))).count_ones() as u64;
        }
        if full_chunks == b8.len() {
            let mask = !1;
            x ^= x << 1;
            x &= mask;
            v += x.count_ones() as u64;
        } else {
            // process x
            let y = u128_from_be_slice(&b8[full_chunks..]);
            v += (x ^ ((x << 1) | (y >> 127))).count_ones() as u64;
            let tail_bits = (b8.len() - full_chunks) * 8;
            let mask = (!0) << (128 - tail_bits + 1);
            let y = y ^ (y << 1);
            v += (y & mask).count_ones() as u64;
        }

        let pi = sample.pop as f64 / n as f64;
        let t = 2.0 * pi * (1.0 - pi);
        let pv =
            erfc(abs(v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)));
        let qv =
            erfc((v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64))) / 2.0;
        TestResult { pv, qv }
    }
}

pub(crate) fn runs_u64(sample: &Sample) -> TestResult {
    let b64 = &sample.b64;
    let n = sample.bit_length;
    let mut v = 1u64;
    for i in 0..b64.len() - 1 {
        let a = b64[i] ^ ((b64[i] << 1) | (b64[i + 1] >> 63));
        v += a.count_ones() as u64;
    }
    let mut last = b64[b64.len() - 1];
    last ^= last << 1;
    let tail_length = n % 64;
    // clear the last bit
    last &= !(1 << ((64 - tail_length) % 64));
    v += last.count_ones() as u64;

    let pi = sample.pop as f64 / n as f64;
    let t = 2.0 * pi * (1.0 - pi);
    let pv = erfc(abs(v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64)));
    let qv = erfc((v as f64 - t * n as f64) / (t * sqrt(2.0 * n as f64))) / 2.0;
    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use super::{super::tests::*, *};
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
        for nbits in 32..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(runs_epsilon(&sample), runs_u64(&sample));
            assert_eq!(runs_epsilon(&sample), runs_u8(&sample));
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
    fn bench_test_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        //  4,928.02 ns/iter
        b.iter(|| {
            test::black_box(runs_u8(&sample));
        });
    }

    #[bench]
    fn bench_test_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 4,313.17 ns/iter
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
