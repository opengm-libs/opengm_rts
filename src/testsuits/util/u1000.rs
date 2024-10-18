use std::{mem::transmute, simd::*};

use num::SimdUint;

use super::u500::*;

// [ff,ff, ff, ff, ...]
// [0, ff, ff, ff, ...]
// [0,  0, ff, ff, ...]
// [..................]
const SHIFT_ELEMENTS_MASK: [u64x16; 16] = {
    let mut res = [u64x16::from_array([!0; 16]); 16];
    let mut buf = [!0; 16];
    let mut i = 0;
    while i < 15 {
        buf[i] = 0;
        res[i + 1] = u64x16::from_array(buf);
        i += 1;
    }
    res
};

#[rustfmt::skip]
const MASK1000: u64x16 = u64x16::from_array(
    {
        let mut a = [!0;16];
        a[15] = (!0) << 24;
        a
    });

#[rustfmt::skip]
const ONE: u64x16 = u64x16::from_array(
    {
        let mut a = [0;16];
        a[0] = 1<<63;
        a
    });

#[derive(Debug, Default, Clone)]
pub struct U1000(u64x16);

impl Bits for U1000 {
    /// from_slice reads 1000 bits from start_pos of b.
    /// start_pos is multiple of 1000.
    fn from_slice(b: &[u8], start_pos: usize) -> Self {
        let b = &b[start_pos / 8..];
        let x:u64x16 = unsafe {
            transmute([
                transmute::<u8x64, u64x8>(u8x64::from_be_bytes(b[64 * 0..64 * 1].try_into().unwrap())).swap_bytes(),
                transmute::<u8x64, u64x8>(u8x64::from_be_bytes(b[64 * 1..64 * 2].try_into().unwrap())).swap_bytes(),
            ])
        };

        // 置零最右12bit
        U1000(x & MASK1000)
    }

    /// from_slice reads 1000 bits from start_pos of b, where b may has
    /// exactly start_pos + 500 bits.
    fn from_slice_exact(b: &[u8], start_pos: usize) -> Self {
        let b = &b[(start_pos - 24) / 8..];
        let x:u64x16 = unsafe {
            transmute([
                transmute::<u8x64, u64x8>(u8x64::from_be_bytes(b[64 * 0..64 * 1].try_into().unwrap())).swap_bytes(),
                transmute::<u8x64, u64x8>(u8x64::from_be_bytes(b[64 * 1..64 * 2].try_into().unwrap())).swap_bytes(),
            ])
        };

        // 右移对齐
        U1000(x).shift_left(24)
    }

    #[inline(always)]
    fn zero() -> Self {
        U1000(u64x16::default())
    }

    // bit 0 = 1
    #[inline(always)]
    fn one() -> Self {
        U1000(ONE)
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
            U1000(a ^ b)
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
            U1000(a ^ b)
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
        self.or(&U1000::one().shift_right(n))
    }
}

impl U1000 {
    // 0 <= n < 64
    #[inline(always)]
    fn shift_elements_left(&self, n: usize) -> Self {
        // let x = self.0 & SHIFT_ELEMENTS_MASK[n];
        // let x = match n {
        //     0 => x,
        //     1 => x.rotate_elements_left::<1>(),
        //     2 => x.rotate_elements_left::<2>(),
        //     3 => x.rotate_elements_left::<3>(),
        //     4 => x.rotate_elements_left::<4>(),
        //     5 => x.rotate_elements_left::<5>(),
        //     6 => x.rotate_elements_left::<6>(),
        //     7 => x.rotate_elements_left::<7>(),
        //     8 => x.rotate_elements_left::<8>(),
        //     9 => x.rotate_elements_left::<9>(),
        //     10 => x.rotate_elements_left::<10>(),
        //     11 => x.rotate_elements_left::<11>(),
        //     12 => x.rotate_elements_left::<12>(),
        //     13 => x.rotate_elements_left::<13>(),
        //     14 => x.rotate_elements_left::<14>(),
        //     15 => x.rotate_elements_left::<15>(),
        //     _ => u64x16::default(),
        // };
        // U1000(x)
        U1000(u64x16_shift_elements_left(&self.0, n))
    }

    // 0 <= n < 16
    #[inline(always)]
    fn shift_elements_right(&self, n: usize) -> Self {
        // let x = self.0;
        // let x = match n {
        //     0 => x,
        //     1 => x.rotate_elements_right::<1>() & SHIFT_ELEMENTS_MASK[1],
        //     2 => x.rotate_elements_right::<2>() & SHIFT_ELEMENTS_MASK[2],
        //     3 => x.rotate_elements_right::<3>() & SHIFT_ELEMENTS_MASK[3],
        //     4 => x.rotate_elements_right::<4>() & SHIFT_ELEMENTS_MASK[4],
        //     5 => x.rotate_elements_right::<5>() & SHIFT_ELEMENTS_MASK[5],
        //     6 => x.rotate_elements_right::<6>() & SHIFT_ELEMENTS_MASK[6],
        //     7 => x.rotate_elements_right::<7>() & SHIFT_ELEMENTS_MASK[7],
        //     8 => x.rotate_elements_right::<8>() & SHIFT_ELEMENTS_MASK[8],
        //     9 => x.rotate_elements_right::<9>() & SHIFT_ELEMENTS_MASK[9],
        //     10 => x.rotate_elements_right::<10>() & SHIFT_ELEMENTS_MASK[10],
        //     11 => x.rotate_elements_right::<11>() & SHIFT_ELEMENTS_MASK[11],
        //     12 => x.rotate_elements_right::<12>() & SHIFT_ELEMENTS_MASK[12],
        //     13 => x.rotate_elements_right::<13>() & SHIFT_ELEMENTS_MASK[13],
        //     14 => x.rotate_elements_right::<14>() & SHIFT_ELEMENTS_MASK[14],
        //     15 => x.rotate_elements_right::<15>() & SHIFT_ELEMENTS_MASK[15],
        //     _ => u64x16::default(),
        // };
        // U1000(x)
        U1000(u64x16_shift_elements_right(&self.0, n))
    }
}


#[inline(always)]
pub(crate) fn u64x16_shift_elements_left(x: &u64x16, n: usize) -> u64x16 {
    let x = x & SHIFT_ELEMENTS_MASK[n%16];
    match n {
        0 => x,
        1 => x.rotate_elements_left::<1>(),
        2 => x.rotate_elements_left::<2>(),
        3 => x.rotate_elements_left::<3>(),
        4 => x.rotate_elements_left::<4>(),
        5 => x.rotate_elements_left::<5>(),
        6 => x.rotate_elements_left::<6>(),
        7 => x.rotate_elements_left::<7>(),
        8 => x.rotate_elements_left::<8>(),
        9 => x.rotate_elements_left::<9>(),
        10 => x.rotate_elements_left::<10>(),
        11 => x.rotate_elements_left::<11>(),
        12 => x.rotate_elements_left::<12>(),
        13 => x.rotate_elements_left::<13>(),
        14 => x.rotate_elements_left::<14>(),
        15 => x.rotate_elements_left::<15>(),
        _ => u64x16::default(),
    }
}

#[inline(always)]
pub(crate) fn u64x16_shift_elements_right(x: &u64x16, n: usize) -> u64x16 {
    match n {
        0 => *x,
        1 => x.rotate_elements_right::<1>() & SHIFT_ELEMENTS_MASK[1],
        2 => x.rotate_elements_right::<2>() & SHIFT_ELEMENTS_MASK[2],
        3 => x.rotate_elements_right::<3>() & SHIFT_ELEMENTS_MASK[3],
        4 => x.rotate_elements_right::<4>() & SHIFT_ELEMENTS_MASK[4],
        5 => x.rotate_elements_right::<5>() & SHIFT_ELEMENTS_MASK[5],
        6 => x.rotate_elements_right::<6>() & SHIFT_ELEMENTS_MASK[6],
        7 => x.rotate_elements_right::<7>() & SHIFT_ELEMENTS_MASK[7],
        8 => x.rotate_elements_right::<8>() & SHIFT_ELEMENTS_MASK[8],
        9 => x.rotate_elements_right::<9>() & SHIFT_ELEMENTS_MASK[9],
        10 => x.rotate_elements_right::<10>() & SHIFT_ELEMENTS_MASK[10],
        11 => x.rotate_elements_right::<11>() & SHIFT_ELEMENTS_MASK[11],
        12 => x.rotate_elements_right::<12>() & SHIFT_ELEMENTS_MASK[12],
        13 => x.rotate_elements_right::<13>() & SHIFT_ELEMENTS_MASK[13],
        14 => x.rotate_elements_right::<14>() & SHIFT_ELEMENTS_MASK[14],
        15 => x.rotate_elements_right::<15>() & SHIFT_ELEMENTS_MASK[15],
        _ => u64x16::default(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mask() {
        for x in SHIFT_ELEMENTS_MASK{
            println!("{:x?}", x);
        }
    }

    #[test]
    fn test_simd() {
        #[rustfmt::skip]
        let b = Vec::from_iter(1..200);
        println!("{:02x?}", b);
        println!();

        let v = U1000::from_slice(&b, 0);
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
    use test::Bencher;

    #[inline(never)]
    fn poptest() {}
    #[bench]
    fn bench_popcount(b: &mut Bencher) {
        let buf = [0xab; 512];
        let u = U1000::from_slice(&buf, 0);
        b.iter(|| {
            for _ in 0..1000000 {
                test::black_box(u.population_parity());
            }
        });
    }
}
