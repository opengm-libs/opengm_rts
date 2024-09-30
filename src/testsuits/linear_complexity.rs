use crate::{igamc, powi, Sample, TestResult};
use crate::testsuits::util::*;

#[cfg(test)]
#[inline(always)]
pub(crate) fn linear_complexity_epsilon(sample: &Sample, m: i32) -> TestResult {
    let mut nu = [0; 7];
    let pi = [
        0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833,
    ];

    let m = m as usize;
    let n = sample.e.len();
    let N = n / m;
    let e = sample.e.as_slice();
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 }; // -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0
        - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
    for i in 0..N {
        let L = berlekamp_massey_epsilon(&e[i * m..(i + 1) * m]) as f64;
        let t = sign * (L - mean) + 2.0 / 9.0;

        if t <= -2.5 {
            nu[0] += 1;
        } else if t > -2.5 && t <= -1.5 {
            nu[1] += 1;
        } else if t > -1.5 && t <= -0.5 {
            nu[2] += 1;
        } else if t > -0.5 && t <= 0.5 {
            nu[3] += 1;
        } else if t > 0.5 && t <= 1.5 {
            nu[4] += 1;
        } else if t > 1.5 && t <= 2.5 {
            nu[5] += 1;
        } else {
            nu[6] += 1;
        }
    }

    let mut chi2 = 0.00;
    for (v, p) in core::iter::zip(nu, pi) {
        let np = N as f64 * p;
        chi2 += powi(v as f64 - np, 2) / np;
    }
    let pv = igamc(3.0, chi2 / 2.0);
    TestResult { pv, qv: pv }
}

//
// set c = c + b*D^e
fn add_shift_epsilon(c: &mut Vec<u8>, b: &Vec<u8>, e: usize) {
    let mut i;
    let mut j;

    if c.len() <= e {
        // |-----c------|
        // |--------e-------|--------b-------|
        i = c.len();
        while i < e {
            c.push(0);
            i += 1;
        }
        for bi in b {
            c.push(*bi);
        }
    } else if c.len() <= e + b.len() {
        // |----------c----------|
        // |--------e-------|--------b-------|
        i = e;
        j = 0;
        while i < c.len() {
            c[i] = (c[i] + b[j]) % 2;
            i += 1;
            j += 1;
        }
        while j < b.len() {
            c.push(b[j]);
            i += 1;
            j += 1;
        }
    } else {
        // |----------------c---------------------|
        // |--------e-------|--------b-------|
        i = e;
        j = 0;
        while j < b.len() {
            c[i] = (c[i] + b[j]) % 2;
            i += 1;
            j += 1;
        }
    }
}

pub(crate) fn berlekamp_massey_epsilon(s: &[u8]) -> usize {
    let n = s.len();
    let mut C: Vec<u8> = Vec::with_capacity(s.len());
    let mut B: Vec<u8> = Vec::with_capacity(s.len());
    let mut T: Vec<u8> = Vec::with_capacity(s.len());

    let mut L = 0;
    let mut m = -1i32;
    C.push(1);
    B.push(1);
    let mut N = 0;
    while N < n {
        let mut d = s[N];
        for i in 1..=L.min(C.len() - 1) {
            d += C[i] * s[N - i];
        }
        d %= 2;

        if d == 1 {
            // T(D) = C(D)
            T.clone_from(&C);
            // C(D) = C(D) + B(D)*D^{N-m}
            let nmm = (N as i32 - m) as usize;
            add_shift_epsilon(&mut C, &B, nmm);

            if L <= N / 2 {
                L = N + 1 - L;
                m = N as i32;
                B.clone_from(&T)
            }
        }
        N += 1;
    }
    L
}

///////////////////////////////////////////////////////////////////////
fn get_bit(s: &[u64], n: usize) -> u8 {
    ((s[n / 64] >> (63 - n % 64)) & 1) as u8
}
fn get_bit_rev(c: &[u64], n: usize) -> u8 {
    let n = c.len() * 64 - n - 1;
    get_bit(c, n)
}

// extract bits from b at start position start_bit, of bit length bit_length.
// zeros may added before the start_bit, returns the total bit length (bit_length + padding zeros).
fn extract(
    b: &[u8],
    start_bit: usize,
    bit_length: usize,
    dst: &mut Vec<u64>,
) -> (usize, usize) {
    dst.clear();
    let leading_zeros = start_bit % 8;
    let start = start_bit / 8;
    let end = (start_bit + bit_length - 1) / 8 + 1;

    let mut i = start;
    while i + 8 < end {
        dst.push(u64::from_be_bytes((&b[i..i + 8]).try_into().unwrap()));
        i += 8;
    }
    if i < end {
        dst.push(u64_from_be_slice(&b[i..end]));
    }

    // clear the leading_zeros bits.
    if leading_zeros > 0 {
        dst[0] &= (1 << (64 - leading_zeros as u32)) - 1;
    }
    (leading_zeros, bit_length)
}
fn extract_borrowed(b64: &[u64],start_bit: usize, bit_length: usize) -> (usize, &[u64]) {
    let leading_zeros = start_bit % 64;
    let start = start_bit / 64;
    let end = (start_bit + bit_length - 1) / 64;

    (leading_zeros, &b64[start..=end])
}

// extract the bits (N-64, N], i.e., bit [N-63, N-62, ..., N].
fn extract_forward_u64(s: &[u64], N: usize) -> u64 {
    // if N is exactly on the boundery.
    if N % 64 == 63 {
        return s[N / 64];
    } else {
        (s[N / 64] >> (63 - N % 64)) | (s[N / 64 - 1] << (N % 64 + 1))
    }
}

/// 线性复杂度检测
pub(crate) fn linear_complexity(sample: &Sample, m: i32) -> TestResult {
    // m = 500, 1000, 5000
    let mut nu = [0; 7];
    let pi = [
        0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833,
    ];

    let m = m as usize;
    let n = sample.bit_length;

    let N = n / m;
    let b64 = sample.b64.as_slice();
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 }; // -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0
        - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
    for i in 0..N {
        let (leadingzeros, S) = extract_borrowed(b64, i*m, m);
        let L = berlekamp_massey(&S, leadingzeros, m) as f64;
        let t = sign * (L - mean) + 2.0 / 9.0;

        if t <= -2.5 {
            nu[0] += 1;
        } else if t > -2.5 && t <= -1.5 {
            nu[1] += 1;
        } else if t > -1.5 && t <= -0.5 {
            nu[2] += 1;
        } else if t > -0.5 && t <= 0.5 {
            nu[3] += 1;
        } else if t > 0.5 && t <= 1.5 {
            nu[4] += 1;
        } else if t > 1.5 && t <= 2.5 {
            nu[5] += 1;
        } else {
            nu[6] += 1;
        }
    }
    let mut chi2 = 0.00;
    for (v, p) in core::iter::zip(nu, pi) {
        let np = N as f64 * p;
        chi2 += powi(v as f64 - np, 2) / np;
    }
    let pv = igamc(3.0, chi2 / 2.0);
    TestResult { pv, qv: pv }
}

//
// set c = c + b*D^e
fn add_shift(c: &mut Vec<u64>, b: &Vec<u64>, e: usize) {
    let left_shift = e % 64;
    match left_shift {
        0 => {
            let mut ci = (c.len() - 1 - e / 64) as isize;
            let mut bi = (b.len() - 1) as isize;
            while ci >= 0 {
                c[ci as usize] ^= b[bi as usize];
                ci -= 1;
                bi -= 1;
            }
        }
        _ => {
            let right_shift = 64 - left_shift;
            let mut ci = (c.len() - 1 - e / 64) as isize;
            let mut bi = (b.len() - 1) as isize;
            let mut last_b = 0;
            while ci >= 0 {
                c[ci as usize] ^= (b[bi as usize] << left_shift) | last_b;
                last_b = b[bi as usize] >> right_shift;
                ci -= 1;
                bi -= 1;
            }
        }
    }
}

// Find L(S) where S start from pos leadingzeros, with nbits bits.
fn berlekamp_massey(S: &[u64], leadingzeros: usize, nbits: usize) -> usize {
    let length = S.len();
    let mut C: Vec<u64> = vec![0; length]; // C is in reverse order, right aligned.
    let mut B: Vec<u64> = vec![0; length];
    let mut T: Vec<u64> = Vec::with_capacity(length);

    let mut L = 0;
    let mut m: i32 = -1;

    // C[D] = 1
    C[length - 1] = 1;
    B[length - 1] = 1;

    let mut N = 0;
    while N < nbits {
        let d = discrepancy(&S, leadingzeros, &C, N, L);
        if d == 1 {
            // T(D) = C(D)
            T.clone_from(&C);

            // C(D) = C(D) + B(D)*D^{N-m}
            let e = (N as i32 - m) as usize;
            add_shift(&mut C, &B, e);

            if L <= N / 2 {
                L = N + 1 - L;
                m = N as i32;
                B.clone_from(&T)
            }
        }
        N += 1;
    }
    L
}

// shift s to buf, then or with c.
fn discrepancy(S: &[u64], s_start: usize, C: &[u64], N: usize, L: usize) -> u8 {
    if false {
        // naive implementation
        let mut d = 0;
        for i in 0..=L {
            let a = get_bit(S, N + s_start - i);
            let b = get_bit_rev(C, i);
            d ^= a & b;
        }
        return d;
    } else {
        let mut pop = 0;
        let nl = s_start + N - L;
        let mut j = C.len() - 1;
        let mut N = s_start + N;

        while N >= nl + 63 {
            let x = extract_forward_u64(S, N);
            pop += (x & C[j]).count_ones();
            N -= 64;
            j -= 1;
        }

        if N >= nl {
            let x;
            if N / 64 != nl / 64 {
                let x1 = S[N / 64] >> (63 - N % 64);
                let x0 = S[nl / 64] & ((1 << (64 - nl % 64)) - 1);
                let x0 = x0 << (N % 64 + 1);
                x = (x0 | x1) & C[j];
            } else {
                x = (S[N / 64] << (nl % 64)) >> (64 - (N - nl + 1));
            }
            pop += (x & C[j]).count_ones();
        }

        (pop & 1) as u8
    }
}

#[cfg(test)]
mod tests {
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_berlekamp_massey() {
        assert_eq!(5, berlekamp_massey_epsilon(&[0, 0, 1, 1, 0, 1, 1, 1, 0]));
        assert_eq!(0, berlekamp_massey_epsilon(&[0, 0, 0]));
        assert_eq!(3, berlekamp_massey_epsilon(&[0, 0, 1]));
        assert_eq!(
            1,
            berlekamp_massey_epsilon(&[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        );
        assert_eq!(2, berlekamp_massey_epsilon(&[1, 1, 0, 1, 1, 0]));
    }
    
    #[test]
    fn test_extract() {
        let b = vec![
            0xf1, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xfe, 0xdc, 0xba,
            0x98, 0x76, 0x54, 0x32, 0x10,
        ];
        let mut dst = Vec::new();
        let n = extract(&b, 2, 65, &mut dst);
        print!("{:?}: ", n);
        for d in &dst {
            print!("{:016x} ", d);
        }
    }

    #[test]
    fn test_extract_forward_u64() {
        let b = vec![0x0123456789abcdef, 0xfedcba9876543210];
        assert_eq!(b[0], extract_forward_u64(&b, 63));
        assert_eq!(0x123456789abcdeff, extract_forward_u64(&b, 67));
    }

    #[test]
    fn test_add_shift() {
        let mut c = vec![0, 0, 1];
        let b = vec![0, 0, 3];
        let e = 127;
        add_shift(&mut c, &b, e);
        println!("{:#?}", c);
    }

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::LinearComplexity);
        let sample: Sample = tv.0.into();
        assert_eq!(linear_complexity_epsilon(&sample, tv.1), tv.2);
    }

    #[test]
    fn test_e() {
        let tv = get_test_vec_e(crate::TestFuncs::LinearComplexity);
        let sample: Sample = E.into();
        assert_eq!(linear_complexity(&sample, tv.0), tv.1);
    }

    #[test]
    fn test_equal() {
        let sample: Sample = E.into();
        for param in [500,1000, 100, 200,301]{
            let res = linear_complexity(&sample, param);
            assert_eq!(res, linear_complexity_epsilon(&sample, param));
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

        // m=500: 15,323,787.50
        // m=1000:20,563,904.20
        b.iter(|| {
            test::black_box(linear_complexity(&sample, 1000));
        });
    }

    #[bench]
    fn bench_test_epsilon(b: &mut Bencher) {
        let sample: Sample = E.into();
        // m=500: 110,799,558.30 ns
        // m=1000:203,625,191.60 ns/iter
        b.iter(|| {
            test::black_box(linear_complexity_epsilon(&sample, 1000));
        });
    }
}
