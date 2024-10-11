use super::{util::*, super::USE_U8};
use crate::{Sample, TestResult};

#[inline(always)]
pub(crate) fn block_frequency(sample: &Sample, m: i32) -> TestResult {
    if USE_U8 {
        block_frequency_u8(sample, m)
    } else {
        block_frequency_u64(sample, m)
    }
}
/// 块内频数检测
#[cfg(test)]
pub(crate) fn block_frequency_epsilon(sample: &Sample, m: i32) -> TestResult {
    let mut v = 0.0;
    let m = m as usize;

    let mut i = 0;
    while i <= sample.e.len() - m {
        let pi = popcount_epsilon(&sample.e[i..i + m]) as f64 / m as f64 - 0.5;
        v += pi * pi;
        i += m;
    }

    let n = (sample.e.len() / m) as f64;
    let pv = igamc(n / 2.0, v * 2.0 * m as f64);
    TestResult { pv, qv: pv }
}

// m = 1000, 10000, 100000
pub(crate) fn block_frequency_u64(sample: &Sample, m: i32) -> TestResult {
    assert!(m >= 8);
    let mut v = 0.0;
    let m = m as usize;
    let N = sample.bit_length / m;
    let mut i = 0;
    while i < N * m {
        let mut pop = 0;
        let start = i / 64;
        let end = (i + m - 1) / 64;
        if end == start {
            pop += clear_lower_bits_u64(
                clear_higher_bits_u64(sample.b64[start], i % 64),
                63 - (i + m - 1) % 64,
            )
            .count_ones() as u64;
        } else {
            pop += clear_higher_bits_u64(sample.b64[start], i % 64).count_ones()
                as u64;
            pop += clear_lower_bits_u64(sample.b64[end], 63 - (i + m - 1) % 64)
                .count_ones() as u64;
            pop += popcount_u64(&sample.b64[start + 1..end]);
        }
        let pi = pop as f64 / m as f64 - 0.5;
        v += pi * pi;
        i += m;
    }

    let n = (sample.bit_length / m) as f64;
    let pv = igamc(n / 2.0, v * 2.0 * m as f64);
    TestResult { pv, qv: pv }
}

// m = 1000, 10000, 100000
pub(crate) fn block_frequency_u8(sample: &Sample, m: i32) -> TestResult {
    assert!(m % 8 == 0);
    let mut v = 0.0;
    let m = m as usize;
    let N = sample.bit_length / m;
    let mut i = 0;
    while i < N * m {
        let mut pop = 0;
        let start = i / 8;
        let end = (i + m) / 8;
        pop += popcount_u8(&sample.b[start..end]);

        let pi = pop as f64 / m as f64 - 0.5;
        v += pi * pi;
        i += m;
    }

    let n = (sample.bit_length / m) as f64;
    let pv = igamc(n / 2.0, v * 2.0 * m as f64);
    TestResult { pv, qv: pv }
}

// use core::arch::aarch64::*;
// use core::simd::*;
// use std::mem::transmute;

// #[inline(always)]
// fn population(x: u8x64) -> u64 {
//     let mut pop = 0;
//     let v: u8x16 = unsafe { transmute(vcntq_u8(transmute(x.resize::<16>(0)))) };
//     pop += v.reduce_sum() as u64;

//     let x = x.rotate_elements_left::<16>();
//     let v: u8x16 = unsafe { transmute(vcntq_u8(transmute(x.resize::<16>(0)))) };
//     pop += v.reduce_sum() as u64;

//     let x = x.rotate_elements_left::<16>();
//     let v: u8x16 = unsafe { transmute(vcntq_u8(transmute(x.resize::<16>(0)))) };
//     pop += v.reduce_sum() as u64;

//     let x = x.rotate_elements_left::<16>();
//     let v: u8x16 = unsafe { transmute(vcntq_u8(transmute(x.resize::<16>(0)))) };
//     pop += v.reduce_sum() as u64;

//     pop
// }

// // m = 1000, 10000, 100000
// pub(crate) fn block_frequency_simd(sample: &Sample, m: i32) -> TestResult {
//     assert!(m >= 8);
//     let mut v = 0.0;
//     let m = m as usize;

//     let N = sample.len() / m;
//     let mut i = 0;
//     let b = &sample.b;

//     let mask = u8x64::from_array([
//         0, 0, 0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
//         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
//         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
//         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
//         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
//         0xff, 0xff, 0xff,
//     ]);

//     while i < N {
//         // load 1000 bits.
//         let v0 = u8x64::from_slice(&b[125 * i..]);
//         let v1 = u8x64::from_slice(&b[125 * i + 64 - 3..]) & mask;

//         let pop =
//             population(v0) + population(v1);

//         let pi = pop as f64 / m as f64 - 0.5;
//         v += pi * pi;
//         i += 1;
//     }

//     let n = (sample.bit_length / m) as f64;
//     let pv = igamc(n / 2.0, v * 2.0 * m as f64);
//     TestResult { pv, qv: pv }
// }

#[cfg(test)]
mod tests {
    use super::{super::tests::*, *};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::BlockFrequency);
        let sample: Sample = tv.0.into();
        assert_eq!(block_frequency_u64(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let _tv = get_test_vec_e(crate::TestFuncs::BlockFrequency);
        let sample: Sample = E.into();
        assert_eq!(
            block_frequency_epsilon(&sample, 1000),
            block_frequency_u8(&sample, 1000)
        );
        assert_eq!(
            block_frequency_u8(&sample, 1000),
            block_frequency_u64(&sample, 1000)
        );
    }

    #[test]
    fn test_equal() {
        let sample: Sample = E.into();
        for m in 10..1000 {
            assert_eq!(
                block_frequency(&sample, m),
                block_frequency_epsilon(&sample, m)
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

    use super::block_frequency_epsilon;

    #[bench]
    fn bench_block_frequency_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 6,719.88 ns/iter
        b.iter(|| {
            test::black_box(block_frequency_u64(&sample, 1000));
        });
    }

    #[bench]
    fn bench_block_frequency_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 6,272.31 ns/iter
        b.iter(|| {
            test::black_box(block_frequency_u8(&sample, 1000));
        });
    }

    #[bench]
    fn bench_block_frequency_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();
        // 12,501.53 ns/iter
        b.iter(|| {
            test::black_box(block_frequency_epsilon(&sample, 10000));
        });
    }
}
