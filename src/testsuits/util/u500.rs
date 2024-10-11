use std::{mem::transmute, simd::*};

use num::SimdUint;

pub trait Bits: Default + Clone{
    /// from_slice reads 500 bits from start_pos of b.
    /// b must have length at least (start_pos + 500) / 8.
    fn from_slice(b: &[u8], start_pos: usize) -> Self;

    /// from_slice reads 500 bits from start_pos of b, where b may has
    /// exactly start_pos + 500 bits.
    fn from_slice_exact(b: &[u8], start_pos: usize) -> Self;

    fn zero() -> Self;

    fn one() -> Self;

    fn population_parity(&self) -> u8;

    fn shift_right(&self, n: usize) -> Self;

    fn shift_left(&self, n: usize) -> Self;

    fn xor(&mut self, rhs: &Self);

    fn or(&mut self, rhs: &Self);

    fn and(&mut self, rhs: &Self);

    fn set_bit(&mut self, n: usize);
}

// [ff,ff, ff, ff, ...]
// [0, ff, ff, ff, ...]
// [0,  0, ff, ff, ...]
// [..................]
const SHIFT_ELEMENTS_MASK: [u64x8; 8] = {
    let mut res = [u64x8::from_array([!0; 8]); 8];
    let mut buf = [!0; 8];
    let mut i = 0;
    while i < 7 {
        buf[i] = 0;
        res[i + 1] = u64x8::from_array(buf);
        i += 1;
    }
    res
};

#[rustfmt::skip]
const MASK500: u64x8 = u64x8::from_array({
    let mut a = [!0;8];
    a[7] = (!0) << 12;
    a
});

#[rustfmt::skip]
const ONE: u64x8 = u64x8::from_array({
    let mut a = [0;8];
    a[0] = 1<<63;
    a
});

#[derive(Debug, Default, Clone)]
pub struct U500(u64x8);

impl Bits for U500 {
    /// from_slice reads 500 bits from start_pos of b.
    /// b must have multiple of 500.
    fn from_slice(b: &[u8], start_pos: usize) -> Self {
        let mut buf = [0; 8];
        let b = &b[start_pos / 8..];
        for i in 0..8 {
            buf[i] =
                u64::from_be_bytes(b[i * 8..(i + 1) * 8].try_into().unwrap());
        }
        let mut x = U500(u64x8::from_array(buf));

        // 左移对齐
        if start_pos % 8 != 0 {
            // 500 | start_pos => start_pos % 8 = 4
            x = x.shift_left(4)
        }

        // 置零最右12bit
        x.0 &= MASK500;
        x
    }

    /// from_slice reads 500 bits from start_pos of b, where b may has
    /// exactly start_pos + 500 bits.
    fn from_slice_exact(b: &[u8], start_pos: usize) -> Self {
        let mut buf = [0; 8];
        let b = &b[(start_pos - 12) / 8..];
        for i in 0..8 {
            buf[i] =
                u64::from_be_bytes(b[i * 8..(i + 1) * 8].try_into().unwrap());
        }
        let x = U500(u64x8::from_array(buf));

        // 右移对齐
        x.shift_left(12)
    }

    #[inline(always)]
    fn zero() -> Self {
        U500(u64x8::default())
    }

    #[inline(always)]
    fn one() -> Self {
        U500(ONE)
    }



    // just returns the parity of population
    #[inline(always)]
    fn population_parity(&self) -> u8 {
        (self.0.reduce_xor().count_ones() & 1) as u8
    }

    #[inline(always)]
    fn shift_right(&self, n: usize) -> Self {
        let x = self.shift_elements_right(n / 64);
        let n = (n % 64) as u64;
        if n > 0 {
            let xx = x.shift_elements_right(1);
            let a = x.0 >> n;
            let b = xx.0 << (64 - n);
            U500(a ^ b)
        } else {
            x
        }
    }

    #[inline(always)]
    fn shift_left(&self, n: usize) -> Self {
        let x = self.shift_elements_left(n / 64);
        let n = (n % 64) as u64;
        if n > 0 {
            let xx = x.shift_elements_left(1);
            let a = x.0 << n;
            let b = xx.0 >> (64 - n);
            U500(a ^ b)
        } else {
            x
        }
    }

    #[inline(always)]
    fn xor(&mut self, rhs: &Self) {
        self.0 ^= rhs.0
    }

    #[inline(always)]
    fn or(&mut self, rhs: &Self) {
        self.0 |= rhs.0
    }

    #[inline(always)]
    fn and(&mut self, rhs: &Self) {
        self.0 &= rhs.0
    }
    #[inline(always)]
    fn set_bit(&mut self, n: usize) {
        self.or(&U500::one().shift_right(n))
    }
}

impl U500 {
    // 0 <= n < 64
    #[inline(always)]
    fn shift_elements_left(&self, n: usize) -> Self {
        let x = self.0 & SHIFT_ELEMENTS_MASK[n];
        let x = match n {
            0 => x,
            1 => x.rotate_elements_left::<1>(),
            2 => x.rotate_elements_left::<2>(),
            3 => x.rotate_elements_left::<3>(),
            4 => x.rotate_elements_left::<4>(),
            5 => x.rotate_elements_left::<5>(),
            6 => x.rotate_elements_left::<6>(),
            7 => x.rotate_elements_left::<7>(),
            _ => u64x8::default(),
        };
        U500(x)
    }

    // 0 <= n < 64
    #[inline(always)]
    fn shift_elements_right(&self, n: usize) -> Self {
        let x = self.0;
        let x = match n {
            0 => x,
            1 => x.rotate_elements_right::<1>() & SHIFT_ELEMENTS_MASK[1],
            2 => x.rotate_elements_right::<2>() & SHIFT_ELEMENTS_MASK[2],
            3 => x.rotate_elements_right::<3>() & SHIFT_ELEMENTS_MASK[3],
            4 => x.rotate_elements_right::<4>() & SHIFT_ELEMENTS_MASK[4],
            5 => x.rotate_elements_right::<5>() & SHIFT_ELEMENTS_MASK[5],
            6 => x.rotate_elements_right::<6>() & SHIFT_ELEMENTS_MASK[6],
            7 => x.rotate_elements_right::<7>() & SHIFT_ELEMENTS_MASK[7],
            _ => u64x8::default(),
        };
        U500(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mask() {
        for i in 0..64 {
            // println!("{:?}", SHIFT_ELEMENTS_MASK[i]);
            println!("{} => self.0.rotate_elements_left::<{}>(),", i, i);
        }
    }

    #[test]
    fn test_simd() {
        #[rustfmt::skip]
        let b = Vec::from_iter(1..100);
        println!("{:02x?}", b);
        println!();

        let v = U500::from_slice(&b, 4);
        println!("{:02x?}\n", v);

        let v = v.shift_elements_right(4);
        println!("{:02x?}\n", v);

        let v = v.shift_elements_left(4);
        println!("{:02x?}\n", v);

        let v = v.shift_right(12);
        println!("{:02x?}\n", v);

        let v = v.shift_left(12);
        println!("{:02x?}\n", v);
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::{test_data::E, Sample};
    use test::Bencher;

    #[bench]
    fn bench_popcount(b: &mut Bencher) {
        let buf = [0xab; 512];
        let u = U500::from_slice(&buf, 0);
        b.iter(|| {
            for _ in 0..1000000 {
                test::black_box(u.population_parity());
            }
        });
    }

    #[bench]
    fn bench_shift(b: &mut Bencher) {
        let buf = [0xab; 512];
        let u = U500::from_slice(&buf, 0);

        // 629,663.28 ns/iter
        b.iter(|| {
            for _ in 0..1000000 {
                test::black_box(u.shift_left(100));
            }
        });
    }
}
