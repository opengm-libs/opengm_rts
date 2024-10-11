use super::{util::*, super::USE_U8};
use crate::{Sample, TestResult};

const L: usize = 7;
const Q: usize = 1280;

/// 通用检测
/// Note that sample must have length greater than 7*1280.
pub(crate) fn universal(sample: &Sample) -> TestResult {
    if USE_U8 {
        universal_u8(sample)
    } else {
        universal_u64(sample)
    }
}

#[cfg(test)]
pub(crate) fn universal_epsilon(sample: &Sample) -> TestResult {
    const E: f64 = 6.1962507;
    const VAR: f64 = 3.125;
    let n = sample.e.len();
    if n < L * Q {
        return TestResult::default();
    }
    let K = n / L - Q;
    let p = 1 << L;
    let mut T = vec![0usize; p];

    let c = 0.7 - 0.8 / (L as f64)
        + (4.0 + 32.0 / (L as f64)) * powf(K as f64, -3.0 / L as f64) / 15.0;
    let sigma = c * sqrt(VAR / K as f64);
    let mut sum = 0.0;
    let mut dec_rep;
    for i in 1..=Q {
        dec_rep = 0;
        for j in 0..L {
            dec_rep += sample.e[(i - 1) * L + j] as usize * (1 << (L - 1 - j));
        }
        T[dec_rep] = i;
    }
    for i in (Q + 1)..=(Q + K) {
        dec_rep = 0;
        for j in 0..L {
            dec_rep += sample.e[(i - 1) * L + j] as usize * (1 << (L - 1 - j));
        }
        sum += ln(i as f64 - (T[dec_rep as usize] as f64)) / ln(2.0);
        T[dec_rep] = i;
    }
    let phi = sum / K as f64;
    let V = (phi - E) / sigma;

    let pv = erfc(abs(V) / SQRT2);
    let qv = erfc(V / SQRT2) / 2.0;
    TestResult { pv, qv }
}

////////////////////////////////////////////////////////////////

// return the number of bits i_i+1_i+2.._i+6
#[inline]
fn get_bits_u64(b64: &[u64], i: usize) -> usize {
    if i % 64 <= 57 {
        ((b64[i / 64] >> 57 - i % 64) & 0x7f) as usize
    } else {
        let h = b64[i / 64] << (7 - (64 - i % 64));
        let l = b64[i / 64 + 1] >> (121 - i % 64);
        ((h | l) & 0x7f) as usize
    }
}

pub(crate) fn universal_u64(sample: &Sample) -> TestResult {
    const E: f64 = 6.1962507;
    const VAR: f64 = 3.125;
    let n = sample.len();
    let b64 = &sample.b64;
    if n < L * Q {
        return TestResult::default();
    }
    let K = n / L - Q;
    let p = 1 << L;
    let mut T = vec![0usize; p];

    let c = 0.7 - 0.8 / (L as f64)
        + (4.0 + 32.0 / (L as f64)) * powf(K as f64, -3.0 / L as f64) / 15.0;
    let sigma = c * sqrt(VAR / K as f64);
    let mut sum = 0.0;
    let mut dec_rep;
    for i in 1..=Q {
        dec_rep = get_bits_u64(b64, (i - 1) * L);
        T[dec_rep] = i;
    }
    for i in (Q + 1)..=(Q + K) {
        dec_rep = get_bits_u64(b64, (i - 1) * L);
        sum += ln(i as f64 - (T[dec_rep as usize] as f64)) / ln(2.0);
        T[dec_rep] = i;
    }
    let phi = sum / K as f64;
    let V = (phi - E) / sigma;

    let pv = erfc(abs(V) / SQRT2);
    let qv = erfc(V / SQRT2) / 2.0;
    TestResult { pv, qv }
}

////////////////////////////////////////////////////////////////

// return the number of bits i i+1 i+2.. i+6
#[inline]
fn get_bits_u8(b8: &[u8], i: usize) -> usize {
    if i % 8 == 0 {
        (b8[i / 8] >> 1) as usize
    } else if i % 8 == 1 {
        (b8[i / 8] & 0x7f) as usize
    } else {
        let h = b8[i / 8] << (i % 8 - 1);
        let l = b8[i / 8 + 1] >> (9 - i % 8);
        ((h | l) & 0x7f) as usize
    }
}

pub(crate) fn universal_u8(sample: &Sample) -> TestResult {
    const E: f64 = 6.1962507;
    const VAR: f64 = 3.125;
    let n = sample.len();
    let b8 = &sample.b;
    if n < L * Q {
        return TestResult::default();
    }
    let K = n / L - Q;
    let p = 1 << L;
    let mut T = vec![0usize; p];

    let c = 0.7 - 0.8 / (L as f64)
        + (4.0 + 32.0 / (L as f64)) * powf(K as f64, -3.0 / L as f64) / 15.0;
    let sigma = c * sqrt(VAR / K as f64);
    let mut sum = 0.0;
    let mut dec_rep;
    for i in 1..=Q {
        dec_rep = get_bits_u8(b8, (i - 1) * L);
        T[dec_rep] = i;
    }
    for i in (Q + 1)..=(Q + K) {
        dec_rep = get_bits_u8(b8, (i - 1) * L);
        sum += ln(i as f64 - (T[dec_rep as usize] as f64)) / ln(2.0);
        T[dec_rep] = i;
    }
    let phi = sum / K as f64;
    let V = (phi - E) / sigma;

    let pv = erfc(abs(V) / SQRT2);
    let qv = erfc(V / SQRT2) / 2.0;
    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use super::{super::tests::*, *};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Universal);
        let sample: Sample = tv.0.into();
        assert_eq!(universal_epsilon(&sample,), tv.2);
        assert_eq!(universal_u64(&sample,), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::Universal);
        let sample: Sample = E.into();
        assert_eq!(tv.1, universal_epsilon(&sample));
        assert_eq!(tv.1, universal_u64(&sample));
        assert_eq!(tv.1, universal_u8(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 9..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(universal_epsilon(&sample), universal_u64(&sample));
            assert_eq!(universal_u8(&sample), universal_u64(&sample));
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
    fn bench_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m= 1,836,168.88 ns/iter
        b.iter(|| {
            test::black_box(universal_epsilon(&sample));
        });
    }

    #[bench]
    fn bench_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=  1,562,443.15 ns/iter
        b.iter(|| {
            test::black_box(universal_u64(&sample));
        });
    }

    #[bench]
    fn bench_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=  1,650,949.58 ns/iter
        b.iter(|| {
            test::black_box(universal_u8(&sample));
        });
    }
}
