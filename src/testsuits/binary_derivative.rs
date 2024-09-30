use super::util::*;
use crate::{Sample, TestResult};

// Note: for k = 2^n - 1, ei = epsilon[i]^epsilon[i+1]^epsilon[i+2]^... ^epsilon[i+k]
//  0: 0        1        2         3          4          5         6  ...
//  1: 01       12       23        34         45         56        67 ...
//  2: 02       13       24        35         46         57        68 ...
//  3: 0123     1234     2345      3456       ...
//  4: 04       15       26        37         ...
//  5: 0145     1256     2367      3478       ...
//  6: 0246     1357     2468      3579       ...
//  7: 01234567 12345678 ...
//  8: 08       19       2(10)     3(11)      ...
//  9: 0189     129(10)
// 10: 028(10)  ...
// 11:
// 12:
// 13:
// 14:
// 15: 0-15     1-16

// K = 3,7,15

/// 二元推导检测
pub(crate) fn binary_derivative(sample: &Sample, k: i32) -> TestResult {
    binary_derivative_u64(sample, k)
    // binary_derivative_epsilon(sample, k)
}

////////////////////////////////////////////////////////////////
#[cfg(test)]
pub(crate) fn binary_derivative_epsilon(sample: &Sample, k: i32) -> TestResult {
    let epsilon = &sample.e;
    let n = epsilon.len();
    let k = k as usize;
    let pv;
    let qv;

    // k is of the form 2^m - 1.
    if 1 == (k + 1).count_ones() {
        // e'[i] = epsilon[i]^epsilon[i+1]^epsilon[i+2]^... ^epsilon[i+7]
        let mut ei = 0u8;
        // epsilon[0] ^ epsilon[1] ^ ... ^ epsilon[k];
        for i in 0..(k + 1) {
            ei ^= epsilon[i];
        }
        let mut sum = ei as i64;
        for i in 1..(n - k) {
            ei = ei ^ epsilon[i - 1] ^ epsilon[i + k];
            sum += ei as i64;
        }
        sum = 2 * sum - (n - k) as i64;
        pv = erfc(abs(sum as f64) / sqrt((n - k) as f64) / SQRT2);
        qv = erfc((sum as f64) / sqrt((n - k) as f64) / SQRT2) / 2.0;
    } else {
        // make a copy
        let mut tmp_epsilon = epsilon.clone();

        for j in 1..(k + 1) {
            for i in 0..(n - j) {
                tmp_epsilon[i] ^= tmp_epsilon[i + 1];
            }
        }
        let mut sum = 0.0;
        for i in 0..(n - k) {
            sum += (2 * (tmp_epsilon[i] as i8) - 1) as f64;
        }

        pv = erfc(abs(sum) / sqrt((n - k) as f64) / SQRT2);
        qv = erfc(sum / sqrt((n - k) as f64) / SQRT2) / 2.0;
    }

    TestResult { pv, qv }
}

////////////////////////////////////////////////////////////////

fn cumulate_binary_derivative<const K: usize>(b64: &[u64], n: usize) -> u64 {
    let mut sum = 0;

    if (n - 1) % 64 >= K {
        // the last b64 has at least K+1 bits.

        for i in 0..b64.len() - 1 {
            let mut x = b64[i];
            for j in 1..=K {
                x ^= (b64[i] << j) | (b64[i + 1] >> (64 - j));
            }

            sum += x.count_ones() as u64;
        }

        // The last u64
        let mut x = b64[b64.len() - 1];
        for j in 1..=K {
            x ^= b64[b64.len() - 1] << j;
        }

        // clear the last K bits
        x >>= (64-n%64)%64 + K;
        sum += x.count_ones() as u64;
        sum
    }else{
        for i in 0..b64.len() - 2 {
            let mut x = b64[i];
            for j in 1..=K {
                x ^= (b64[i] << j) | (b64[i + 1] >> (64 - j));
            }

            sum += x.count_ones() as u64;
        }

        // b64[-2] || b64[-1]
        let mut x = b64[b64.len() - 2];
        for j in 1..=K {
            x ^= (b64[b64.len() - 2] << j) | (b64[b64.len() - 1] >> (64 - j));
        }

        // clear the last K - n%64 bits
        x >>= K - n%64;
        sum += x.count_ones() as u64;
        sum
    }
}
pub(crate) fn binary_derivative_u64(sample: &Sample, k: i32) -> TestResult {
    assert!(k == 3 || k == 7 || k == 15);

    let b64 = &sample.b64;
    let n = sample.len();
    let k = k as usize;
    let pv;
    let qv;

    
    // k is of the form 2^m - 1.
    // e'[i] = epsilon[i]^epsilon[i+1]^epsilon[i+2]^... ^epsilon[i+7]
    let sum = match k {
            3 => {
                cumulate_binary_derivative::<3>(b64, n)
            }
            7 => {
                cumulate_binary_derivative::<7>(b64, n)
            }
            15 => {
                cumulate_binary_derivative::<15>(b64, n)
            }
            _ => unreachable!(),
    };

    let sum = 2 * sum as i64 - (n - k) as i64;
    pv = erfc(abs(sum as f64) / sqrt((n - k) as f64) / SQRT2);
    qv = erfc((sum as f64) / sqrt((n - k) as f64) / SQRT2) / 2.0;

    TestResult { pv, qv }
}

#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::BinaryDerivative);
        let sample: Sample = tv.0.into();
        assert_eq!(binary_derivative_epsilon(&sample, tv.1), tv.2);
        assert_eq!(binary_derivative_u64(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::BinaryDerivative);
        let sample: Sample = E.into();
        assert_eq!(tv.1, binary_derivative_epsilon(&sample, tv.0));
        assert_eq!(tv.1, binary_derivative_u64(&sample, tv.0));
    }

    #[test]
    fn test_equal() {
        for nbits in 16..2000 {
            let sample: Sample = E[..nbits * 8].into();
            for k in [3, 7, 15] {
                assert_eq!(
                    binary_derivative_epsilon(&sample, k),
                    binary_derivative_u64(&sample, k)
                );
            }
        }
    }
}


#[cfg(test)]
mod bench {
    extern crate test;
    use crate::test_data::E;

    use super::*;
    use super::{Sample};
    use test::Bencher;

    const K:i32 = 3;

    #[bench]
    fn bench_test(b: &mut Bencher) {
        let sample: Sample = E.into();

        // K = 3: 820,142.70 ns/iter
        // K = 7: 823,033.86 ns/iter
        // K = 15: 825,439.59 ns/iter
        b.iter(|| {
            test::black_box(binary_derivative_epsilon(&sample,K));
        });
    }

    #[bench]
    fn bench_test_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();

        // K = 3: 15,437.97 ns/iter
        // K = 7: 35,166.07 ns/iter 
        // K = 15: 75,817.23 ns/iter
        b.iter(|| {
            test::black_box(binary_derivative_u64(&sample, K));
        });
    }
}
