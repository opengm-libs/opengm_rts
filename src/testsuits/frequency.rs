use super::util::*;
use crate::{Sample, TestResult};

/// 频度检测
pub(crate) fn frequency(sample: &Sample) -> TestResult {
    let s = sample.pop as i64 * 2 - (sample.len() as i64);
    let pv = erfc(abs(s as f64) / sqrt((sample.len() * 2) as f64));
    let qv = erfc(s as f64 / sqrt((sample.len() * 2) as f64)) / 2.0;
    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use crate::{test_data::E, Sample};
    use super::{*, super::tests::*};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Frequency);
        let sample: Sample = tv.0.into();
        assert_eq!(frequency(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::Frequency);
        let sample: Sample = E.into();
        assert_eq!(tv.1, frequency(&sample));
    }
}
