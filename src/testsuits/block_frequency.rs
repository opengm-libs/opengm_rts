use super::util::*;
use crate::{Sample, TestResult};

#[inline(always)]
pub(crate) fn block_frequency(sample: &Sample, m: i32) -> TestResult {
    block_frequency_u64(sample, m)
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
pub(crate) fn block_frequency_u8(sample: &Sample, m: i32) -> TestResult {
    assert!(m >= 8);

    if m % 8 != 0 {
        // if true{
        // also works for m % 8 == 0.
        let mut v = 0.0;
        let m = m as usize;
        let mut i = 0;
        while i <= sample.b.len() * 8 - m {
            let mut pop = 0u64;
            if i % 8 != 0 {
                pop += (sample.b[i / 8] & lower_bits_mask(8 - i % 8))
                    .count_ones() as u64;
            };
            if (i + m) % 8 != 0 {
                pop += (sample.b[(i + m) / 8] & higher_bits_mask((i + m) % 8))
                    .count_ones() as u64;
            }
            pop += popcount(&sample.b[(i + 7) / 8..(i + m) / 8]) as u64;

            let pi = pop as f64 / m as f64 - 0.5;
            v += pi * pi;
            i += m;
        }

        let n = (sample.b.len() * 8 / m) as f64;
        let pv = igamc(n / 2.0, v * 2.0 * m as f64);
        TestResult { pv, qv: pv }
    } else {
        let mut v = 0.0;
        let m = m as usize;
        let m8 = m >> 3;

        let mut i = 0;
        while i <= sample.b.len() - m8 {
            let pi = popcount(&sample.b[i..i + m8]) as f64 / m as f64 - 0.5;
            v += pi * pi;
            i += m8;
        }

        let n = (sample.b.len() * 8 / m) as f64;
        let pv = igamc(n / 2.0, v * 2.0 * m as f64);
        TestResult { pv, qv: pv }
    }
}

// m = 1000, 10000, 100000
pub(crate) fn block_frequency_u64(sample: &Sample, m: i32) -> TestResult {
    assert!(m >= 8);
    let mut v = 0.0;
    let m = m as usize;
    let N = sample.bit_length/m;
    let mut i = 0;
    while i < N * m {
        let mut pop = 0;
        let start = i/64;
        let end = (i+m-1)/64;
        if end == start{
            pop += clear_lower_bits_u64(clear_higher_bits_u64(sample.b64[start], i%64),63-(i+m-1)%64).count_ones() as u64;
        }else{
            pop += clear_higher_bits_u64(sample.b64[start], i%64).count_ones() as u64;
            pop += clear_lower_bits_u64(sample.b64[end], 63-(i+m-1)%64).count_ones()as u64;
            pop += popcount_u64(&sample.b64[start+1..end]);
        }
        let pi = pop as f64 / m as f64 - 0.5;
        v += pi * pi;
        i += m;
    }

    let n = (sample.bit_length / m) as f64;
    let pv = igamc(n / 2.0, v * 2.0 * m as f64);
    TestResult { pv, qv: pv }
}

#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::BlockFrequency);
        let sample: Sample = tv.0.into();
        assert_eq!(block_frequency_u64(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::BlockFrequency);
        let sample: Sample = E.into();
        assert_eq!(block_frequency(&sample, tv.0), tv.1);
    }

    #[test]
    fn test_equal() {
        let sample: Sample = E.into();
        for m in 10..1000{
            assert_eq!(block_frequency(&sample, m), block_frequency_epsilon(&sample, m));
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
    fn bench_block_frequency_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 2,960.34
        b.iter(|| {
            test::black_box(block_frequency_u8(&sample, 10000));
        });
    }

    #[bench]
    fn bench_block_frequency_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 3,484.13 ns/iter
        b.iter(|| {
            test::black_box(block_frequency_u64(&sample, 10000));
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
