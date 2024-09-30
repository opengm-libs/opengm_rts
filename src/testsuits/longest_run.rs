use super::util::*;
use crate::{get_bit_unchecked, Sample, TestResult};

/// 块内最大0游程检测
pub(crate) fn longest_run0(sample: &Sample) -> TestResult {
    longest_run0_u64(sample)
}

/// 块内最大0游程检测
pub(crate) fn longest_run1(sample: &Sample) -> TestResult {
    longest_run1_u64(sample)
}


// ref: B.7 of GM/T 0005-2021
// Assume n >= 128
#[inline(always)]
fn get_params(n: usize) -> (i32, i32, [i32; 7], [f64; 7]) {
    let K;
    let M;
    let mut V = [0; 7];
    let mut pi = [0.0; 7];

    // ref: B.7 of GM/T 0005-2021
    if n < 6272 {
        K = 3;
        M = 8;
        V[0] = 1;
        V[1] = 2;
        V[2] = 3;
        V[3] = 4;
        pi[0] = 0.2148;
        pi[1] = 0.3672;
        pi[2] = 0.2305;
        pi[3] = 0.1875;
    } else if n < 750000 {
        K = 5;
        M = 128;
        V[0] = 4;
        V[1] = 5;
        V[2] = 6;
        V[3] = 7;
        V[4] = 8;
        V[5] = 9;
        pi[0] = 0.1174;
        pi[1] = 0.2430;
        pi[2] = 0.2494;
        pi[3] = 0.1752;
        pi[4] = 0.1027;
        pi[5] = 0.1124;
    } else {
        K = 6;
        M = 10000;
        V[0] = 10;
        V[1] = 11;
        V[2] = 12;
        V[3] = 13;
        V[4] = 14;
        V[5] = 15;
        V[6] = 16;
        pi[0] = 0.086632;
        pi[1] = 0.208201;
        pi[2] = 0.248419;
        pi[3] = 0.193913;
        pi[4] = 0.121458;
        pi[5] = 0.068011;
        pi[6] = 0.073366;
    }
    (K, M, V, pi)
}

////////////////////////////////////////////////////////////////

#[inline(always)]
fn longest_run_internal_epsilon<const BIT: u8>(e: &[u8], K: i32, M: i32, V: &[i32], pi: &[f64]) -> f64 {
    let n = e.len();
    let N = n / (M as usize);
    let mut nu = [0; 7]; // K at most 6

    for i in 0..N {
        let block = &e[i * (M as usize)..];
        let mut longestRun = 0;
        let mut currentRun = 0;
        for j in 0..(M as usize) {
            if block[j] == BIT {
                currentRun += 1;
                if longestRun < currentRun {
                    longestRun = currentRun;
                }
            } else {
                currentRun = 0;
            }
        }
        nu[(saturating(longestRun, V[0], V[K as usize]) - V[0]) as usize] += 1;
    }

    let mut chi2 = 0.0;
    for i in 0..(K + 1) {
        let N = N as f64;
        let nui = nu[i as usize] as f64;
        let pii = pi[i as usize];
        chi2 += (nui - N * pii) * (nui - N * pii) / (N * pii);
    }
    igamc(K as f64 / 2.0, chi2 / 2.0)
}


/// 块内最大0游程检测
#[cfg(test)]
pub(crate) fn longest_run0_epsilon(sample: &Sample) -> TestResult {
    if sample.e.len() < 128 {
        panic!("longest run test: n too short, at least 128\n");
    }

    let (K, M, V, pi) = get_params(sample.e.len());

    let pv = longest_run_internal_epsilon::<0>(&sample.e, K, M, &V, &pi);

    TestResult {
        pv,
        qv: pv,
    }
}

/// 块内最大0游程检测
#[cfg(test)]
pub(crate) fn longest_run1_epsilon(sample: &Sample) -> TestResult {
    if sample.e.len() < 128 {
        panic!("longest run test: n too short, at least 128\n");
    }

    let (K, M, V, pi) = get_params(sample.e.len());

    let pv = longest_run_internal_epsilon::<1>(&sample.e, K, M, &V, &pi);

    TestResult {
        pv,
        qv: pv,
    }
}

////////////////////////////////////////////////////////////////////////



#[inline(always)]
fn longest_run_internal_u64<const BIT: u8>(b64: &[u64], n:usize, K: i32, M: i32, V: &[i32], pi: &[f64]) -> f64 {
    let N = n / (M as usize);
    let mut nu = [0; 7]; // K at most 6

    for i in 0..N {
        let block_start = i * (M as usize);
        let mut longestRun = 0;
        let mut currentRun = 0;
        for j in 0..(M as usize) {
            if get_bit_unchecked(b64, block_start + j) == BIT {
                currentRun += 1;
                if longestRun < currentRun {
                    longestRun = currentRun;
                }
            } else {
                currentRun = 0;
            }
        }
        nu[(saturating(longestRun, V[0], V[K as usize]) - V[0]) as usize] += 1;
    }

    let mut chi2 = 0.0;
    for i in 0..(K + 1) {
        let N = N as f64;
        let nui = nu[i as usize] as f64;
        let pii = pi[i as usize];
        chi2 += (nui - N * pii) * (nui - N * pii) / (N * pii);
    }
    igamc(K as f64 / 2.0, chi2 / 2.0)
}

/// 块内最大0游程检测
pub(crate) fn longest_run0_u64(sample: &Sample) -> TestResult {
    if sample.len() < 128 {
        panic!("longest run test: n too short, at least 128\n");
    }

    let (K, M, V, pi) = get_params(sample.len());

    let pv = longest_run_internal_u64::<0>(&sample.b64, sample.bit_length, K, M, &V, &pi);

    TestResult {
        pv,
        qv: pv,
    }
}

/// 块内最大0游程检测
pub(crate) fn longest_run1_u64(sample: &Sample) -> TestResult {
    if sample.len() < 128 {
        panic!("longest run test: n too short, at least 128\n");
    }

    let (K, M, V, pi) = get_params(sample.len());

    let pv = longest_run_internal_u64::<1>(&sample.b64, sample.bit_length, K, M, &V, &pi);

    TestResult {
        pv,
        qv: pv,
    }
}


#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::LongestRun0);
        let sample: Sample = tv.0.into();
        assert_eq!(longest_run0_epsilon(&sample), tv.2);
        assert_eq!(longest_run0_u64(&sample), tv.2);
        
        let tv = get_test_vec(crate::TestFuncs::LongestRun1);
        let sample: Sample = tv.0.into();
        assert_eq!(longest_run1_epsilon(&sample), tv.2);
        assert_eq!(longest_run1_u64(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::LongestRun0);
        let sample: Sample = E.into();
        assert_eq!(tv.1, longest_run0_epsilon(&sample));
        assert_eq!(tv.1, longest_run0_u64(&sample));

        let tv = get_test_vec_e(crate::TestFuncs::LongestRun1);
        let sample: Sample = E.into();
        assert_eq!(tv.1, longest_run1_epsilon(&sample));
        assert_eq!(tv.1, longest_run1_u64(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 128/8..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(longest_run0_epsilon(&sample), longest_run0_u64(&sample));
            assert_eq!(longest_run1_epsilon(&sample), longest_run1_u64(&sample));
        }
    }
}
