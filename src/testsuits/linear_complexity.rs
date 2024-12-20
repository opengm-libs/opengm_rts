use std::cmp::max;
use std::mem::{self, swap};

// use u1000::U1000;
// use u500::U500;

use crate::testsuits::util::*;
use crate::{igamc, powi, Sample, TestResult, USE_U8};

pub(crate) fn linear_complexity(sample: &Sample, m: i32) -> TestResult {
    if USE_U8 {
        // when compute m = 1000, we are also computed the first half 500
        // bits's complexity. So we can save half the computation for m = 500
        // when sample.len() = 1 million.
        if m == 500 || m == 1000{
            linear_complexity5001000_u8(sample, m)
        }else{
            linear_complexity_u8(sample, m)
        }
    } else {
        linear_complexity_u64(sample, m)
    }

    // if m == 500 {
    //     linear_complexity_simd::<U500, 500>(sample)
    // } else if m == 1000 {
    //     linear_complexity_simd::<U1000, 1000>(sample)
    // } else {
    //     linear_complexity_u8(sample, m)
    // }
}

#[cfg(test)]
#[inline(always)]
pub(crate) fn linear_complexity_epsilon(sample: &Sample, m: i32) -> TestResult {
    let mut nu = [0; 7];
    let pi = [0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833];

    let m = m as usize;
    let n = sample.e.len();
    let N = n / m;
    let e = sample.e.as_slice();
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 }; // -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0 - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
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
            // d += C[i] * s[N - i];
            d ^= C[i] & s[N - i];
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
fn extract(b: &[u8], start_bit: usize, bit_length: usize, dst: &mut Vec<u64>) -> (usize, usize) {
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

fn extract_borrowed_u64(b64: &[u64], start_bit: usize, bit_length: usize) -> (usize, &[u64]) {
    let leading_zeros = start_bit % 64;
    let start = start_bit / 64;
    let end = (start_bit + bit_length - 1) / 64;

    (leading_zeros, &b64[start..=end])
}

fn extract_borrowed_u8(out: &mut [u64], b8: &[u8], start_bit: usize, bit_length: usize) -> usize {
    let leading_zeros = start_bit % 64;
    let mut start = start_bit / 64 * 8;
    out[0] = u64::from_be_bytes(b8[start..start + 8].try_into().unwrap());
    start += 8;
    let mut i = 1;
    let mut count = 64 - leading_zeros; // have read count bits.
    while count < bit_length {
        out[i] = u64::from_be_bytes(b8[start..start + 8].try_into().unwrap());
        i += 1;
        start += 8;
        count += 64;
    }

    leading_zeros
}

// extract the bits (N-64, N], i.e., bit [N-63, N-62, ..., N].
#[inline(always)]
fn extract_forward_u64(s: &[u64], N: usize) -> u64 {
    if true {
        (s[N / 64] >> (63 - N % 64)) | ((s[(N - 63) / 64] << 1) << (N % 64))
    } else {
        // if N is exactly on the boundery.
        if N % 64 == 63 {
            return s[N / 64];
        } else {
            (s[N / 64] >> (63 - N % 64)) | (s[N / 64 - 1] << (N % 64 + 1))
        }
    }
}

/// 线性复杂度检测
pub(crate) fn linear_complexity_u8(sample: &Sample, m: i32) -> TestResult {
    // m = 500, 1000, 5000
    let mut nu = [0; 7];
    let pi = [0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833];

    let m = m as usize;
    let n = sample.bit_length;

    let N = n / m;
    let b8 = &sample.b;
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 }; // -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0 - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
    let mut S = vec![0; (m + 63) / 64 + 1];
    for i in 0..N {
        // let (leadingzeros, S) = extract_borrowed_u64(b64, i * m, m);
        let leadingzeros = extract_borrowed_u8(&mut S, b8, i * m, m);
        let L = berlekamp_massey_u64(&S, leadingzeros, m) as f64;
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


/// 线性复杂度检测
pub(crate) fn linear_complexity5001000_u8(sample: &Sample, m: i32) -> TestResult {
    // m = 500, 1000, 5000
    let mut nu = [0; 7];
    let pi = [0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833];

    let m = m as usize;
    let n = sample.bit_length;

    let N = n / m;
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 }; // -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0 - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
    let complexity = get_complexity(sample, m);
    for L in complexity {
        let L = L as f64;
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

/// 线性复杂度检测
pub(crate) fn linear_complexity_u64(sample: &Sample, m: i32) -> TestResult {
    // m = 500, 1000, 5000
    let mut nu = [0; 7];
    let pi = [0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833];

    let m = m as usize;
    let n = sample.bit_length;

    let N = n / m;
    let b64 = sample.b64.as_slice();
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 }; // -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0 - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
    for i in 0..N {
        let (leadingzeros, S) = extract_borrowed_u64(b64, i * m, m);
        let L = berlekamp_massey_u64(&S, leadingzeros, m) as f64;
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

// set c = c + b*D^e
fn add_shift_u64(c: &mut Vec<u64>, b: &Vec<u64>, e: usize) {
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

//
// set c = c + b*D^e
#[inline(always)]
fn add_shift_u64_with_length(c: &mut [u64], b: &[u64], e: usize, clen: usize, blen: usize) {
    let left_shift = e % 64;
    let mut op_bits = max(clen, blen + e) as isize;
    match left_shift {
        0 => {
            let mut ci = (c.len() - 1 - e / 64) as isize;
            let mut bi = (b.len() - 1) as isize;
            while (ci >= 0) && (op_bits > 0) {
                c[ci as usize] ^= b[bi as usize];
                ci -= 1;
                bi -= 1;
                op_bits -= 64;
            }
        }
        _ => {
            let right_shift = 64 - left_shift;
            let mut ci = (c.len() - 1 - e / 64) as isize;
            let mut bi = (b.len() - 1) as isize;
            let mut last_b = 0;
            while (ci >= 0) && (op_bits > 0) {
                c[ci as usize] ^= (b[bi as usize] << left_shift) | last_b;
                last_b = b[bi as usize] >> right_shift;
                ci -= 1;
                bi -= 1;
                op_bits -= 64;
            }
        }
    }
}

// Find L(S) where S start from pos leadingzeros, with nbits bits.
fn berlekamp_massey_u64(S: &[u64], leadingzeros: usize, nbits: usize) -> usize {
    let length = S.len();
    let mut C: Vec<u64> = vec![0; length]; // C is in reverse order, right aligned.
    let mut B: Vec<u64> = vec![0; length];
    let mut T: Vec<u64> = vec![0; length];

    let mut L = 0;
    let mut m: i32 = -1;

    // C[D] = 1
    C[length - 1] = 1;
    B[length - 1] = 1;
    let mut blen = 1;

    let mut N = 0;
    while N < nbits {
        let d = discrepancy_u64(&S, leadingzeros, &C, N, L);
        if d == 1 {
            // T(D) = C(D)
            T[length - (L + 64) / 64..].copy_from_slice(&C[length - (L + 64) / 64..]);
            // tlen = L+1;

            // C(D) = C(D) + B(D)*D^{N-m}
            let e = (N as i32 - m) as usize;
            // add_shift_u64(&mut C, &B, e);
            // deg(C) = L
            add_shift_u64_with_length(&mut C, &B, e, L + 1, blen);

            if L <= N / 2 {
                // B.copy_from_slice(&T);
                swap(&mut B, &mut T);
                blen = L + 1;
                m = N as i32;
                L = N + 1 - L;
            }
        }
        N += 1;
    }
    L
}

////////////////////////////////////////////////////////////////////////////////////////////////

// only for m = 500 or 1000
fn get_complexity(sample: &Sample, m: usize) -> Vec<u16> {
    if let Ok(mut complexities) = sample.complexities.lock() {
        if let Some(v) = &mut complexities[m / 500 - 1] {
            return mem::take(v);
        }
    }

    let (complexity500, complexity1000) = compute_complexity_1000(sample);

    if let Ok(mut complexities) = sample.complexities.lock() {
        if m == 500 {
            complexities[1] = Some(complexity1000);
            return complexity500;
        }else{
            complexities[0] = Some(complexity500);
            return complexity1000;
        }
    }else{
        unreachable!()
    }

}

// compute the complexities for m = 500, 1000
fn compute_complexity_1000(sample: &Sample) -> (Vec<u16>, Vec<u16>) {
    let n = sample.len();

    // assume n mod 1000 < 500, that is N500 = 2N1000.
    let N500 = n / 500;
    let N1000 = n / 1000;
    debug_assert!(N500 == 2 * N1000);

    let mut complexity500 = Vec::with_capacity(N500);
    let mut complexity1000 = Vec::with_capacity(N1000);

    let b8 = &sample.b;
    let mut S = vec![0; 17]; // 16*64 = 1024
    for i in 0..N1000 {
        let leadingzeros = extract_borrowed_u8(&mut S, b8, i * 1000, 1000);
        let (l1, l2, l3) = berlekamp_massey1000_u64(&S, leadingzeros, 1000);
        complexity500.push(l1 as u16);
        complexity500.push(l2 as u16);
        complexity1000.push(l3 as u16);
    }
    return (complexity500, complexity1000);
}

// Find L(S) where S start from pos leadingzeros, with nbits bits.
// the special case nbits = 1000, returns
// the first L(S[..500]) and L(S)
fn berlekamp_massey2_u64(S: &[u64], leadingzeros: usize, nbits: usize) -> (usize, usize) {
    let length = S.len();
    let mut C: Vec<u64> = vec![0; length]; // C is in reverse order, right aligned.
    let mut B: Vec<u64> = vec![0; length];
    let mut T: Vec<u64> = vec![0; length];

    let mut L = 0;
    let mut m: i32 = -1;

    // C[D] = 1
    C[length - 1] = 1;
    B[length - 1] = 1;
    let mut blen = 1;

    let mut N = 0;
    let mut res1 = 0;
    while N < nbits {
        if N == nbits / 2 {
            res1 = L;
        }
        let d = discrepancy_u64(&S, leadingzeros, &C, N, L);
        if d == 1 {
            // T(D) = C(D)
            T[length - (L + 64) / 64..].copy_from_slice(&C[length - (L + 64) / 64..]);
            // tlen = L+1;

            // C(D) = C(D) + B(D)*D^{N-m}
            let e = (N as i32 - m) as usize;
            // add_shift_u64(&mut C, &B, e);
            // deg(C) = L
            add_shift_u64_with_length(&mut C, &B, e, L + 1, blen);

            if L <= N / 2 {
                // B.copy_from_slice(&T);
                swap(&mut B, &mut T);
                blen = L + 1;
                m = N as i32;
                L = N + 1 - L;
            }
        }
        N += 1;
    }
    (res1, L)
}

// Find L(S) where S start from pos leadingzeros, with nbits bits.
// the special case nbits = 1000, returns L(S[..500]) L(S[500..]) and L(S)
fn berlekamp_massey1000_u64(S: &[u64], leadingzeros: usize, nbits: usize) -> (usize, usize, usize) {
    let length = S.len();
    let mut C: Vec<u64> = vec![0; length]; // C is in reverse order, right aligned.
    let mut B: Vec<u64> = vec![0; length];
    let mut T: Vec<u64> = vec![0; length];

    let mut L = 0;
    let mut m: i32 = -1;

    // C[D] = 1
    C[length - 1] = 1;
    B[length - 1] = 1;
    let mut blen = 1;

    let mut N = 0;
    let mut res1 = 0;
    while N < nbits {
        if N == nbits / 2 {
            res1 = L;
        }
        let d = discrepancy_u64(&S, leadingzeros, &C, N, L);
        if d == 1 {
            // T(D) = C(D)
            T[length - (L + 64) / 64..].copy_from_slice(&C[length - (L + 64) / 64..]);
            // tlen = L+1;

            // C(D) = C(D) + B(D)*D^{N-m}
            let e = (N as i32 - m) as usize;
            // add_shift_u64(&mut C, &B, e);
            // deg(C) = L
            add_shift_u64_with_length(&mut C, &B, e, L + 1, blen);

            if L <= N / 2 {
                // B.copy_from_slice(&T);
                swap(&mut B, &mut T);
                blen = L + 1;
                m = N as i32;
                L = N + 1 - L;
            }
        }
        N += 1;
    }
    let res2 = berlekamp_massey_u64(&S[(leadingzeros + 500)/64..], (leadingzeros + 500) % 64, 500);
    (res1, res2, L)
}

// d = C_0S_N + C_1S_{N-1} + ... + C_LS_{N-L} % 2
#[inline(always)]
fn discrepancy_u64(S: &[u64], s_start: usize, C: &[u64], N: usize, L: usize) -> u8 {
    if false {
        // naive implementation
        let mut d = 0;
        for i in 0..=L {
            let a = get_bit(S, N + s_start - i);
            let b = get_bit_rev(C, i);
            d ^= a & b;
        }
        return d;
    } else if false {
        let mut pop = 0;
        let nl = s_start + N - L;
        let mut j = C.len() - 1;
        let mut N = s_start + N;

        while N >= nl + 63 {
            let x = extract_forward_u64(S, N);
            pop ^= (x & C[j]).count_ones();
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
            pop ^= (x & C[j]).count_ones();
        }

        (pop & 1) as u8
    } else {
        let mut pop = 0;
        let N_minus_L = s_start + N - L;
        let N = s_start + N;

        let mut c_idx = C.len() - 1;
        // S[N-L..N] exactly on the boundary
        let rshift = 63 - N % 64;
        if rshift == 0 {
            let mut s_idx = (N / 64) as i64;
            while s_idx >= (N_minus_L / 64) as i64 {
                pop ^= C[c_idx] & S[s_idx as usize];
                c_idx -= 1;
                s_idx -= 1;
            }
        } else {
            let mut s_idx = (N / 64) as i64;
            let lshift = 64 - rshift;
            while s_idx > (N_minus_L / 64) as i64 {
                pop ^= C[c_idx] & ((S[s_idx as usize] >> rshift) | (S[s_idx as usize - 1] << lshift));
                s_idx -= 1;
                c_idx -= 1;
            }
            pop ^= C[c_idx] & (S[s_idx as usize] >> rshift);
        }

        (pop.count_ones() & 1) as u8
    }
}

#[cfg(test)]
mod tests {
    use super::{super::tests::*, *};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_berlekamp_massey() {
        assert_eq!(5, berlekamp_massey_epsilon(&[0, 0, 1, 1, 0, 1, 1, 1, 0]));
        assert_eq!(0, berlekamp_massey_epsilon(&[0, 0, 0]));
        assert_eq!(3, berlekamp_massey_epsilon(&[0, 0, 1]));
        assert_eq!(1, berlekamp_massey_epsilon(&[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]));
        assert_eq!(2, berlekamp_massey_epsilon(&[1, 1, 0, 1, 1, 0]));
    }

    #[test]
    fn test_extract() {
        let b = vec![
            0xf1, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10,
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
        add_shift_u64(&mut c, &b, e);
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
        // let tv = get_test_vec_e(crate::TestFuncs::LinearComplexity);
        let sample: Sample = E.into();
        // for m in [500, 1000, 5000] {
        for m in [500, 1000] {
            println!("{}", m);
            // FIXME: panic when m = 5000
            let c = linear_complexity_epsilon(&sample, m);
            let a = linear_complexity(&sample, m);
            // let b = linear_complexity_u64(&sample, m);

            // assert_eq!(a,b);
            assert_eq!(a, c);
        }
    }

    #[test]
    fn test_equal() {
        let sample: Sample = E.into();
        for param in [500, 1000, 100, 200, 301] {
            let res = linear_complexity_u8(&sample, param);
            assert_eq!(res, linear_complexity_epsilon(&sample, param));
        }
    }

    #[test]
    fn test_berlekamp_massey2() {
        let sample: Sample = E.into();
        let mut S = vec![0; (1000 + 63) / 64 + 1];
        let start = 100;
        let leading_zeros = extract_borrowed_u8(&mut S, &sample.b, start, 1000);
        let (l1, l2) = berlekamp_massey2_u64(&S, leading_zeros, 1000);
        let leading_zeros = extract_borrowed_u8(&mut S, &sample.b, start, 500);
        let l3 = berlekamp_massey_u64(&S, leading_zeros, 500);
        assert_eq!(l1, l3);
        println!("{} {} {}", l1, l2, l3);
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::{test_data::E, Sample};
    use test::Bencher;

    #[bench]
    fn bench_linear_complexity_u64(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=500: 15,102,658.30
        // m=1000:20,563,904.20
        // m=5000:61,863,154.20
        b.iter(|| {
            test::black_box(linear_complexity_u64(&sample, 5000));
        });
    }

    #[bench]
    fn bench_linear_complexity_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // m=500:  13,285,025.00
        // m=1000: 15,692,979.10
        // m=5000: 33,255,095.90
        b.iter(|| {
            // test::black_box(linear_complexity_u8(&sample, 500));
            // test::black_box(linear_complexity_u8(&sample, 1000));
            test::black_box(linear_complexity_u8(&sample, 1000));
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
