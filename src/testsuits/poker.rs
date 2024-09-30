use super::util::*;
use crate::{Sample, TestResult};

const MAX_M: usize = 8;

/// 扑克检测, only for m = 4 or 8.

#[inline(always)]
pub(crate) fn poker(sample: &Sample, m: i32) -> TestResult {
    poker_u64(sample, m)
}


/// 扑克检测
#[cfg(test)]
pub(crate) fn poker_epsilon(sample: &Sample, m: i32) -> TestResult {
    assert!(m == 4 || m == 8);

    let m = m as usize;
    let n = sample.e.len() / m;
    let power: usize = 1 << m;
    let mut ni = vec![0; power];

    for i in 0..n {
        let mut idx = 0_usize;
        for j in 0..m {
            idx = 2 * idx + sample.e[i * m + j] as usize;
        }
        ni[idx] += 1;
    }

    let sum = ni.iter().map(|x| x * x).sum::<usize>();
    let v = (power as f64) / (n as f64) * (sum as f64) - (n as f64);
    let pv = igamc((power as f64 - 1.0) / 2.0, v / 2.0);
    TestResult { pv, qv: pv }
}

/// 扑克检测
pub(crate) fn poker_u8(sample: &Sample, m: i32) -> TestResult {
    assert!(m == 4 || m == 8);
    // We assume sample always has bits of multiples of 8.
    assert!(sample.bit_length %8==0);

    let m = m as usize;
    let n = sample.bit_length / m;
    let power: usize = 1 << m;
    let mut ni = vec![0; power];
    if m == 4 {
        for b in &sample.b {
            ni[(*b >> 4) as usize] += 1;
            ni[(*b & 0x0f) as usize] += 1;
        }
    } else {
        for b in &sample.b {
            ni[*b as usize] += 1;
        }
    }

    let sum = ni.iter().map(|x| x * x).sum::<usize>();
    let v = (power as f64) / (n as f64) * (sum as f64) - (n as f64);
    let pv = igamc((power as f64 - 1.0) / 2.0, v / 2.0);
    TestResult { pv, qv: pv }
}


/// 扑克检测
pub(crate) fn poker_u64(sample: &Sample, m: i32) -> TestResult {
    assert!(m == 4 || m == 8);

    let m = m as usize;
    let n = sample.bit_length / m;
    let power: usize = 1 << m;
    let mut ni = vec![0; power];

    if sample.bit_length % 64 == 0{
        if m == 4 {
            for b in &sample.b64 {
                let mut x = *b;
                for _ in 0..16{
                    ni[(x & 0x0f) as usize] += 1;
                    x >>= 4;
                }
            }
        } else {
            for b in &sample.b64 {
                let mut x = *b;
                for _ in 0..8{
                    ni[x as u8 as usize] += 1;
                    x >>= 8;
                }
            }
        }
    }else{
        let b64_length = sample.b64.len();
        if m == 4 {
            for b in &sample.b64[..b64_length-1] {
                let mut x = *b;
                for _ in 0..16{
                    ni[(x & 0x0f) as usize] += 1;
                    x >>= 4;
                }
            }
        } else {
            for b in &sample.b64[..b64_length-1] {
                let mut x = *b;
                for _ in 0..8{
                    ni[x as u8 as usize] += 1;
                    x >>= 8;
                }
            }
        }
        // handle the last one
        let mut tail_length = sample.bit_length % 64;
        let mut x = sample.b64[b64_length-1];
        while tail_length >= m{
            ni[(x>>(64-m)) as usize] += 1;
            x <<= m;
            tail_length -= m;
        }
    }


    let sum = ni.iter().map(|x| x * x).sum::<usize>();
    let v = (power as f64) / (n as f64) * (sum as f64) - (n as f64);
    let pv = igamc((power as f64 - 1.0) / 2.0, v / 2.0);
    TestResult { pv, qv: pv }
}


#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Poker);
        let sample: Sample = tv.0.into();
        assert_eq!(poker(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::Poker);
        let sample: Sample = E.into();
        assert_eq!(tv.1, poker(&sample, tv.0));
    }

    #[test]
    fn test_equal() {
        let sample: Sample = E[..9992].into();
        assert_eq!(poker_epsilon(&sample, 4), poker_u8(&sample, 4));
        assert_eq!(poker_epsilon(&sample, 8), poker_u8(&sample, 8));
        assert_eq!(poker_epsilon(&sample, 4), poker_u64(&sample, 4));
        assert_eq!(poker_epsilon(&sample, 8), poker_u64(&sample, 8));
    }
}

#[cfg(test)]
mod bench{
    extern crate test;
    use test::Bencher;
    use super::*;
    use crate::{test_data::E, Sample};
    #[bench]
    fn bench_poker_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=4: 101,742.19 ns
        // m=8: 50,177.50 ns
        b.iter(|| {
            test::black_box(poker_u8(&sample, 8));
        });
    }

    #[bench]
    fn bench_poker_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=4: 84,803.82 ns/iter
        // m=8: 37,354.41 ns/iter
        b.iter(|| {
            test::black_box(poker_u64(&sample, 4));
        });
    }

    #[bench]
    fn bench_poker_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();
        // m=4: 550,129.68 ns
        // m=8: 432,989.58 ns
        b.iter(|| {
            test::black_box(poker_epsilon(&sample, 8));
        });
    }

}