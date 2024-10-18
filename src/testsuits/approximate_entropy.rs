use super::{util::*, super::USE_U8};
use crate::{Sample, TestResult};
// m = 2, 5, 7
const MAX_M: usize = 8;

/// 近似熵检测
pub(crate) fn approximate_entropy(sample: &Sample, m: i32) -> TestResult {
    if USE_U8 {
        approximate_entropy_u8(sample, m)
    } else {
        approximate_entropy_u64(sample, m)
    }
}

////////////////////////////////////////////////////////////////////////

#[cfg(test)]
pub(crate) fn approximate_entropy_epsilon(
    sample: &Sample,
    m: i32,
) -> TestResult {
    let n = sample.e.len();
    let apen = phi_epsilon(sample.e.as_slice(), m)
        - phi_epsilon(sample.e.as_slice(), m + 1);

    let v = 2.0 * (n as f64) * (ln(2.0) - apen);
    let pv = igamc(powi(2.0, m - 1), v / 2.0);
    TestResult { pv, qv: pv }
}

fn phi_epsilon(e: &[u8], m: i32) -> f64 {
    if m == 0 {
        return 0.0;
    }
    if m as usize > MAX_M {
        panic!("Invalid m:{}", m)
    }

    let n = e.len();
    let m = m as usize;
    let mut probs = [0usize; 1 << MAX_M];

    // buf = {e[n-m+1], ..., e[n-1], e[0], e[1], e[m-2]}
    let mut buf = [0u8; 2 * MAX_M]; // buf.len = m-1 + m-1, only 2m-2 are need.
    for i in 0..(2 * m - 2) {
        buf[i] = e[(n - m + 1 + i) % n];
    }

    for w in e.windows(m) {
        probs[get_partten_epsilon(w)] += 1;
    }

    for w in buf[..(2 * m - 2)].windows(m) {
        probs[get_partten_epsilon(w)] += 1;
    }

    let mut sum = 0.0;
    for i in 0..(1 << m) {
        if probs[i] == 0 {
            continue;
        }
        let c = probs[i] as f64 / n as f64;
        sum += c * ln(c);
    }
    sum
}

fn get_partten_epsilon(e: &[u8]) -> usize {
    let mut k = 0;
    for i in e {
        k <<= 1;
        k |= *i as usize;
    }
    k
}

////////////////////////////////////////////////////////////////////////

pub(crate) fn approximate_entropy_u64(sample: &Sample, m: i32) -> TestResult {
    assert!(m >= 0 && m as usize <= MAX_M);

    let n = sample.len();
    let phi_m = phi_u64(&sample, m);
    let phi_m1 = phi_u64(&sample, m + 1);
    let apen = phi_m - phi_m1;

    let v = 2.0 * (n as f64) * (ln(2.0) - apen);
    let pv = igamc(powi(2.0, m - 1), v / 2.0);
    TestResult { pv, qv: pv }
}

// m = 2, 5
fn phi_u64(sample: &Sample, m: i32) -> f64 {
    let patterns = overlapping_patterns_u64(sample, m as usize);
    let mut sum = 0.0;
    for x in patterns.iter() {
        if *x == 0 {
            continue;
        }
        let c = *x as f64 / sample.bit_length as f64;
        sum += c * ln(c);
    }
    sum
}

////////////////////////////////////////////////////////////////////////

pub(crate) fn approximate_entropy_u8(sample: &Sample, m: i32) -> TestResult {
    assert!(m >= 0 && m as usize <= MAX_M);

    let n = sample.len();
    let phi_m = phi_u8(&sample, m);
    let phi_m1 = phi_u8(&sample, m + 1);
    let apen = phi_m - phi_m1;

    let v = 2.0 * (n as f64) * (ln(2.0) - apen);
    let pv = igamc(powi(2.0, m - 1), v / 2.0);
    TestResult { pv, qv: pv }
}

// m = 2, 5
fn phi_u8(sample: &Sample, m: i32) -> f64 {
    let patterns = overlapping_patterns_u8(sample, m as usize);
    let mut sum = 0.0;
    for x in patterns.iter() {
        if *x == 0 {
            continue;
        }
        let c = *x as f64 / sample.bit_length as f64;
        sum += c * ln(c);
    }
    sum
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
        let tv = get_test_vec(crate::TestFuncs::ApproximateEntropy);
        let sample: Sample = tv.0.into();
        assert_eq!(approximate_entropy_epsilon(&sample, tv.1), tv.2);
        assert_eq!(approximate_entropy_u64(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::ApproximateEntropy);
        let sample: Sample = E.into();
        assert_eq!(tv.1, approximate_entropy_epsilon(&sample, tv.0));
        assert_eq!(tv.1, approximate_entropy_u64(&sample, tv.0));
        assert_eq!(tv.1, approximate_entropy_u8(&sample, tv.0));
    }

    #[test]
    fn test_equal() {
        for nbits in 128 / 8..1000 {
            let sample: Sample = E[..nbits * 8].into();
            for d in [2, 5, 7] {
                assert_eq!(
                    approximate_entropy_epsilon(&sample, d),
                    approximate_entropy_u64(&sample, d)
                );

                assert_eq!(
                    approximate_entropy_u8(&sample, d),
                    approximate_entropy_u64(&sample, d)
                );
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
    const M: i32 = 5;
    #[bench]
    fn bench_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m= 7,615,637.35 ns/iter
        b.iter(|| {
            test::black_box(approximate_entropy_epsilon(&sample, M));
        });
    }

    #[bench]
    fn bench_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=  1,418.13 ns/iter
        b.iter(|| {
            test::black_box(approximate_entropy_u64(&sample, M));
        });
    }

    #[bench]
    fn bench_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=  1,436.72 ns/iter
        b.iter(|| {
            test::black_box(approximate_entropy_u8(&sample, M));
        });
    }
}
