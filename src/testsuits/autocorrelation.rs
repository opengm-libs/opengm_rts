use super::{super::USE_U8, util::*};
use crate::{Sample, TestResult};

/// 自相关检测
pub(crate) fn autocorrelation(sample: &Sample, d: i32) -> TestResult {
    if USE_U8 {
        autocorrelation_u8(sample, d)
    } else {
        autocorrelation_u64(sample, d)
    }
}

////////////////////////////////////////////////////////////////
#[cfg(test)]
pub(crate) fn autocorrelation_epsilon(sample: &Sample, d: i32) -> TestResult {
    let epsilon = &sample.e;
    let n = epsilon.len();
    let d = d as usize;

    let mut sum = 0i64;
    for i in 0..(n - d) {
        sum += (epsilon[i] ^ epsilon[i + d]) as i64;
    }

    let v = (2.0 * sum as f64 - (n - d) as f64) / sqrt((n - d) as f64);
    let pv = erfc(abs(v) / SQRT2);
    let qv = erfc(v / SQRT2) / 2.0;

    TestResult { pv, qv }
}

////////////////////////////////////////////////////////////////

pub(crate) fn autocorrelation_u64(sample: &Sample, d: i32) -> TestResult {
    assert!(d == 1 || d == 2 || d == 4 || d == 8 || d == 16 || d == 32);

    let b64 = &sample.b64;
    let n = sample.len();
    let d = d as usize;

    let mut sum = 0;

    if (n - 1) % 64 >= d {
        for i in 0..b64.len() - 1 {
            let x = b64[i] ^ ((b64[i] << d) | (b64[i + 1] >> (64 - d)));
            sum += x.count_ones() as u64;
        }

        let mut x = b64[b64.len() - 1];
        x ^= x << d;
        // clear the tail bits
        x >>= (64 - (n % 64)) % 64 + d;
        sum += x.count_ones() as u64;
    } else {
        for i in 0..b64.len() - 2 {
            let x = b64[i] ^ ((b64[i] << d) | (b64[i + 1] >> (64 - d)));
            sum += x.count_ones() as u64;
        }

        let mut x = b64[b64.len() - 2] ^ ((b64[b64.len() - 2] << d) | (b64[b64.len() - 1] >> (64 - d)));
        // clear the tail bits
        x >>= d - (n % 64);
        sum += x.count_ones() as u64;
    }

    let v = (2.0 * sum as f64 - (n - d) as f64) / sqrt((n - d) as f64);
    let pv = erfc(abs(v) / SQRT2);
    let qv = erfc(v / SQRT2) / 2.0;

    TestResult { pv, qv }
}

////////////////////////////////////////////////////////////////

pub(crate) fn autocorrelation_u8(sample: &Sample, d: i32) -> TestResult {
    assert!(d == 1 || d == 2 || d == 4 || d == 8 || d == 16 || d == 32);
    let b8 = &sample.b;
    let n = sample.len();
    let d = d as usize;
    assert!(b8.len() >= 8);

    let mut sum = 0;

    // iterator bench will more faster: 4,969.89 ns/iter
    let mut iter = sample.iter_u64();
    let mut x = iter.next().unwrap();
    for y in iter {
        // x || y
        let z = x ^ ((x << d) | (y >> (64 - d)));
        sum += z.count_ones() as u64;
        x = y;
    }
    let (tail, tail_bits) = sample.tail_u64();

    // 5,608.38 ns/iter 
    // let full_chunks = b8.len() & (!7);
    // let mut x = u64::from_be_bytes((&b8[..8]).try_into().unwrap());
    // for chunk in b8[8..full_chunks].chunks_exact(8) {
    //     let y = u64::from_be_bytes(chunk.try_into().unwrap());

    //     // x || y
    //     let z = x ^ ((x << d) | (y >> (64 - d)));
    //     sum += z.count_ones() as u64;

    //     x = y;
    // }

    // // leave x || tail to go
    // let tail = if full_chunks < b8.len() {
    //     u64_from_be_slice(&b8[full_chunks..])
    // } else {
    //     0
    // };
    // // 0 <= tail_bits < 64
    // let tail_bits = (b8.len() - full_chunks) * 8;

    // 64-K <= bits_to_go < 128 - K
    let bits_to_go = 64 + tail_bits - d;
    if bits_to_go <= 64 {
        // only x need to process
        let mut z = x ^ ((x << d) | (tail >> (64 - d)));
        z = clear_lower_bits_u64(z, 64 - bits_to_go);
        sum += z.count_ones() as u64;
    } else {
        // process x
        let mut z = x ^ ((x << d) | (tail >> (64 - d)));
        sum += z.count_ones() as u64;

        // process tail
        z = tail ^ (tail << d);
        z = clear_lower_bits_u64(z, 128 - bits_to_go);
        sum += z.count_ones() as u64;
    }

    let v = (2.0 * sum as f64 - (n - d) as f64) / sqrt((n - d) as f64);
    let pv = erfc(abs(v) / SQRT2);
    let qv = erfc(v / SQRT2) / 2.0;

    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        test_data::E,
        testsuits::tests::{get_test_vec, get_test_vec_e},
        Sample,
    };

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Autocorrelation);
        let sample: Sample = tv.0.into();
        assert_eq!(autocorrelation_epsilon(&sample, tv.1), tv.2);
        assert_eq!(autocorrelation_u64(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::Autocorrelation);
        let sample: Sample = E.into();
        assert_eq!(tv.1, autocorrelation_epsilon(&sample, tv.0));
        assert_eq!(tv.1, autocorrelation_u64(&sample, tv.0));
        assert_eq!(tv.1, autocorrelation_u8(&sample, tv.0));
    }

    #[test]
    fn test_equal() {
        for nbits in 128 / 8..1000 {
            let sample: Sample = E[..nbits * 8].into();
            for d in [1, 2, 4, 8, 16, 32] {
                assert_eq!(autocorrelation_epsilon(&sample, d), autocorrelation_u64(&sample, d));

                assert_eq!(autocorrelation_u8(&sample, d), autocorrelation_u64(&sample, d));
            }
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
    fn bench_u64(b: &mut Bencher) {
        let sample: Sample = E.into();
        b.iter(|| {
            test::black_box(autocorrelation_u64(&sample, 16));
        });
    }

    #[bench]
    fn bench_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=16: 4,952.14 ns/iter
        b.iter(|| {
            test::black_box(autocorrelation_u8(&sample, 16));
        });
    }

    #[bench]
    fn bench_test_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();
        b.iter(|| {
            test::black_box(autocorrelation_epsilon(&sample, 16));
        });
    }
}
