use super::util::*;
use crate::Sample;

pub(crate) fn linear_complexity(sample: &Sample, m: i32) -> f64 {
    let mut nu = [0; 7];
    let pi = [
        0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833,
    ];

    let m = m as usize;
    let n = sample.e.len();
    let N = n / m;
    let e = sample.e.as_slice();
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 };// -1^m
    let mean = m as f64 / 2.0 + (9.0 - sign) / 36.0
        - 1.0 / powi(2.0, m as i32) * (m as f64 / 3.0 + 2.0 / 9.0);
    for i in 0..N {
        let L = berlekamp_massey(&e[i * m..(i + 1) * m]) as f64;
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
    for (v,p) in core::iter::zip(nu, pi) {
        let np = N as f64 * p;
        chi2 += powi(v as f64 - np, 2) / np;
    }
    igamc(3.0, chi2 / 2.0)
}

//
// set c = c + b*D^e
fn add_shift(c: &mut Vec<u8>, b: &Vec<u8>, e: usize) {
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

fn berlekamp_massey(s: &[u8]) -> usize {
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
            add_shift(&mut C, &B, nmm);

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_berlekamp_massey() {
        assert_eq!(5, berlekamp_massey(&[0, 0, 1, 1, 0, 1, 1, 1, 0]));
        assert_eq!(0, berlekamp_massey(&[0, 0, 0]));
        assert_eq!(3, berlekamp_massey(&[0, 0, 1]));
        assert_eq!(1, berlekamp_massey(&[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]));
        assert_eq!(2, berlekamp_massey(&[1, 1, 0, 1, 1, 0]));
    }
}