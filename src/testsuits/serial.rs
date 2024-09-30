use super::util::*;
use crate::{Sample, TestResult};

const MAX_M: usize = 7;
const MIN_M: usize = 2;

/// 重叠子序列p_value1检测
#[inline(always)]
pub(crate) fn serial1(sample: &Sample, m: i32) -> TestResult {
    serial1_u64(sample, m)
}

/// 重叠子序列p_value2检测
#[inline(always)]
pub(crate) fn serial2(sample: &Sample, m: i32) -> TestResult {
    serial2_u64(sample, m)
}

////////////////////////////////////////////////////////////////

fn psi2_epsilon(e: &[u8], m: i32) -> f64 {
    let n = e.len();
    if m <= 0 {
        return 0.0;
    }

    let length = 1 << m;
    let mask = length - 1;
    let mut bucket = [0u64; 1 << MAX_M];

    let mut k = 0;
    let mut i: usize = 0;
    while i < m as usize {
        k = (k << 1) | (e[i] as usize);
        i += 1;
    }
    bucket[k] += 1;

    for _ in 1..n {
        k = (k << 1) | (e[i] as usize);
        k &= mask;
        bucket[k] += 1;
        i += 1;
        if i == n {
            i = 0;
        }
    }
    // Note that a uint sum may overflow.
    let sum: f64 = bucket[..1 << m].iter().map(|x| ((*x) * (*x)) as f64).sum();
    sum / (n as f64) * (1 << m) as f64 - n as f64
}

/// 重叠子序列p_value1检测
#[cfg(test)]
pub(crate) fn serial1_epsilon(sample: &Sample, m: i32) -> TestResult {
    assert!(m <= MAX_M as i32 && m >= MIN_M as i32);

    let p0 = psi2_epsilon(&sample.e, m);
    let p1 = psi2_epsilon(&sample.e, m - 1);
    let del1 = p0 - p1;
    let pv = igamc(powi(2.0, m - 2), del1 / 2.0);

    TestResult { pv, qv: pv }
}

/// 重叠子序列p_value2检测
#[cfg(test)]
pub(crate) fn serial2_epsilon(sample: &Sample, m: i32) -> TestResult {
    assert!(m <= MAX_M as i32 && m >= MIN_M as i32);

    let p0 = psi2_epsilon(&sample.e, m);
    let p1 = psi2_epsilon(&sample.e, m - 1);
    let p2 = psi2_epsilon(&sample.e, m - 2);
    let del2 = p0 - 2.0 * p1 + p2;

    let pv = igamc(powi(2.0, m - 3), del2 / 2.0);

    TestResult { pv, qv: pv }
}

////////////////////////////////////////////////////////////////////////

// m = 3,5,7 in GM/T 0005-2021.
pub(crate) fn psi2_u64(sample: &Sample, M: usize) -> f64 {
    assert!(sample.len() >= 2);
    if M == 0{
        return 0.0
    }

    let n = sample.len();
    let patterns = overlapping_patterns(sample, M);

    let sum: f64 = patterns.iter().map(|x| ((*x) * (*x)) as f64).sum();
    sum / (n as f64) * ((1 << M) as f64) - (n as f64)
}

/// 重叠子序列p_value1检测
pub(crate) fn serial1_u64(sample: &Sample, m: i32) -> TestResult {
    // m = 3,5,7

    assert!(m <= MAX_M as i32 && m >= MIN_M as i32);

    let p0 = psi2_u64(sample, m as usize);
    let p1 = psi2_u64(sample, m as usize - 1);
    let del1 = p0 - p1;
    let pv = igamc(powi(2.0, m - 2), del1 / 2.0);

    TestResult { pv, qv: pv }
}

/// 重叠子序列p_value2检测
pub(crate) fn serial2_u64(sample: &Sample, m: i32) -> TestResult {
    // note 2 <= m <= 10
    assert!(m <= MAX_M as i32 && m >= MIN_M as i32);

    let p0 = psi2_u64(sample, m as usize);
    let p1 = psi2_u64(sample, m as usize - 1);
    let p2 = psi2_u64(sample, m as usize - 2);

    let del2 = p0 - 2.0 * p1 + p2;
    let pv = igamc(powi(2.0, m - 3), del2 / 2.0);
    TestResult { pv, qv: pv }
}

#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Serial1);
        let sample: Sample = tv.0.into();
        assert_eq!(serial1_u64(&sample, tv.1), tv.2);
        assert_eq!(serial1_epsilon(&sample, tv.1), tv.2);

        let tv = get_test_vec(crate::TestFuncs::Serial2);
        let sample: Sample = tv.0.into();
        assert_eq!(serial2_u64(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv1 = get_test_vec_e(crate::TestFuncs::Serial1);
        let tv2 = get_test_vec_e(crate::TestFuncs::Serial2);
        let sample: Sample = E.into();

        assert_eq!(tv1.1, serial1_epsilon(&sample, tv1.0));
        assert_eq!(tv1.1, serial1_u64(&sample, tv1.0));
        assert_eq!(tv2.1, serial2_u64(&sample, tv2.0));
    }

    #[test]
    fn test_equal() {
        for nbits in 9..1200 {
            // the sample's bits length must be multiple of 8.
            let sample: Sample = E[..nbits * 8].into();
            for m in MIN_M..=MAX_M {
                let m = m as i32;
                assert_eq!(
                    serial1_epsilon(&sample, m),
                    serial1_u64(&sample, m)
                );
                assert_eq!(
                    serial2_epsilon(&sample, m),
                    serial2_u64(&sample, m)
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

    #[bench]
    fn bench_test(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=3: 2,232,047.90 ns/iter
        // m=5: 1,148,277.05 ns/iter
        // m=7: 886,577.10 ns/iter
        b.iter(|| {
            test::black_box(serial1_u64(&sample, 7));
        });
    }

    #[bench]
    fn bench_test_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();
        // m=3: 2,562,231.20 ns/iter
        // m=5: 2,416,622.90 ns/iter
        // m=7: 2,438,622.90 ns/iter
        b.iter(|| {
            test::black_box(serial1_epsilon(&sample, 7));
        });
    }
}
