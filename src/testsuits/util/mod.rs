mod math;
pub(crate) use math::*;

use crate::Sample;

// returns the sum of epsilon[0] + ... + epsilon[n-1]
// assumes that epsilon[i] = 0 or 1.
pub(crate) fn popcount_epsilon(e: &[u8]) -> u64 {
    let mut sum = 0u64;

    // 2040 = 255*8
    for chunk in e.chunks_exact(255 * 8) {
        let tmp_sum: u64 = chunk
            .chunks_exact(8)
            .into_iter()
            .map(|x| u64::from_ne_bytes(x.try_into().unwrap()))
            .sum();
        sum += (0..8).map(|i| (tmp_sum >> 8 * i) & 0xff).sum::<u64>();
    }

    let remainder = &e[e.len() - e.len() % (255 * 8)..];
    let tmp_sum: u64 = remainder
        .chunks_exact(8)
        .into_iter()
        .map(|x| u64::from_ne_bytes(x.try_into().unwrap()))
        .sum();
    sum += (0..8).map(|i| (tmp_sum >> 8 * i) & 0xff).sum::<u64>();

    let remainder = &remainder[remainder.len() - remainder.len() % 8..];
    sum += remainder.iter().sum::<u8>() as u64;

    sum
}

// each bytes in b contains 8 bits
pub(crate) fn popcount(b: &[u8]) -> u64 {
    let mut pop = 0;

    let full_blocks16 = b.len() & (!15);
    for chunk in b[..full_blocks16].chunks_exact(16) {
        pop += u128::from_ne_bytes(chunk.try_into().unwrap()).count_ones();
    }

    let full_blocks8 = b.len() & (!7);
    for chunk in b[full_blocks16..full_blocks8].chunks_exact(8) {
        pop += u64::from_ne_bytes(chunk.try_into().unwrap()).count_ones();
    }

    let mut tail = 0;
    for i in b[full_blocks8..].iter() {
        tail <<= 8;
        tail |= (*i) as u64;
    }

    pop += tail.count_ones();
    pop as u64
}

// each bytes in b contains 8 bits
pub(crate) fn popcount_u64(b64: &[u64]) -> u64 {
    let mut pop = 0;

    for b in b64 {
        pop += (*b).count_ones();
    }

    pop as u64
}

#[inline(always)]
pub(crate) fn saturating<T: core::cmp::PartialOrd>(n: T, min: T, max: T) -> T {
    if n >= max {
        max
    } else if n <= min {
        min
    } else {
        n
    }
}

#[inline(always)]
pub(crate) fn saturating_ceil<T: core::cmp::PartialOrd>(n: T, max: T) -> T {
    if n >= max {
        max
    }else{
        n
    }
}

// return 0x000...111 of n "1" .
#[inline(always)]
pub(crate) fn lower_bits_mask(n: usize) -> u8 {
    (!0) >> (8 - n)
}

// return 0x11...0 of n "1" .
#[inline(always)]
pub(crate) fn higher_bits_mask(n: usize) -> u8 {
    (!0) << (8 - n)
}

// return 0x000...111 of n "1" .
// 0 <= n <= 64
#[inline(always)]
pub(crate) fn lower_bits_mask_u64(n: usize) -> u64 {
    (!0) >> (64 - n)
}

// return 0x11...0 of n "1" .
#[inline(always)]
pub(crate) fn higher_bits_mask_u64(n: usize) -> u64 {
    (!0) << (64 - n)
}

// set the higher n bits to zero.
// 0 <= n < 64
#[inline(always)]
pub(crate) fn clear_higher_bits_u64(x: u64, n: usize) -> u64 {
    (x << n) >> n
}

// set the lower n bits to zero.
// 0 <= n < 64
#[inline(always)]
pub(crate) fn clear_lower_bits_u64(x: u64, n: usize) -> u64 {
    (x >> n) << n
}

// return 0x0..1..1...0 with leadingzeros '0', and length of '1'.
#[inline(always)]
pub(crate) fn bits_mask_u64(leadingzeros: usize, ones_length: usize) -> u64 {
    ((1 << ones_length) - 1) << (64 - leadingzeros - ones_length)
}

// construct a u64 from a slice, where d.len() <= 8, padding zeros at LSB.
// For example, d = [0x12,0x34,0x56,0x78], then returns 0x1234567800000000.
#[inline(always)]
pub(crate) fn u64_from_be_slice(d: &[u8]) -> u64 {
    let mut x = 0;
    let mut i = 0;
    while i < d.len() {
        x = (x << 8) | d[i] as u64;
        i += 1;
    }
    x <<= (8 - i) * 8;
    x
}

// 处理一个u64,统计重叠子序列模式, tail为上一个u64未处理的最后m-1比特.
#[inline(always)]
fn overlapping_patterns_process_u64(
    bucket: &mut [u64],
    x: u64,
    mask: u64,
    tail: u64,
    tail_mask: u64,
    m: u64,
) -> u64 {
    // process the last m bits of x
    let new_tail = x & tail_mask;
    let mut x = x;
    for _ in 0..(m - 1) {
        bucket[(x & mask) as usize] += 1;
        x >>= 1;
    }

    x |= tail << (64 - (m - 1));

    for _ in 0..(64 - (m - 1)) {
        bucket[(x & mask) as usize] += 1;
        x >>= 1;
    }
    new_tail
}

// 重叠m子序列计数
// 重叠子序列和近似熵测试使用.
// 重叠子序列m的可能取值为(m,m-1,m-2): 0,1,2,3,4,5,6,7
// 近似熵m可能取值:2,5,7
// 参见:GM/T 0005-2021, 附录A
pub(crate) fn overlapping_patterns(sample: &Sample, M: usize) -> Vec<u64> {
    if M == 0 {
        return Vec::new();
    }

    if M == 1 {
        return vec![sample.len() as u64 - sample.pop, sample.pop];
    }

    if let Ok(pat) = sample.patterns.lock() {
        if pat.len() >= M + 1 {
            if let Some(v) = &pat[M] {
                return v.clone();
            }
        }
    }

    let nbits = sample.len();
    let b64 = &sample.b64;
    let mask = (1u64 << M) - 1;
    let tail_mask = (1u64 << (M - 1)) - 1;
    let mut bucket = vec![0u64; 1 << M];

    // The first
    let mut x = b64[0];
    let mut tail = x & tail_mask;
    for _ in (M - 1)..64 {
        bucket[(x & mask) as usize] += 1;
        x >>= 1;
    }

    let last_b64_bits = nbits % 64;
    if last_b64_bits == 0 {
        // The hole b64 are epsilon
        for x in b64[1..].iter() {
            tail = overlapping_patterns_process_u64(
                &mut bucket,
                *x,
                mask,
                tail,
                tail_mask,
                M as u64,
            );
        }

        // last m-1 bits tail
        let mut x = (tail << (M - 1)) | b64[0] >> (64 - (M - 1));
        for _ in 0..M - 1 {
            bucket[(x & mask) as usize] += 1;
            x >>= 1;
        }
    } else {
        // the last one of b64 is not full.
        for x in b64[1..b64.len() - 1].iter() {
            tail = overlapping_patterns_process_u64(
                &mut bucket,
                *x,
                mask,
                tail,
                tail_mask,
                M as u64,
            );
        }

        // tail || b64[-1][63..63-last_b64_bits] || b64[0][..m-1]
        let mut x: u128 = (b64[0] >> (64 - (M - 1))) as u128;
        x |= ((b64[b64.len() - 1] >> (64 - last_b64_bits)) << (M - 1)) as u128;
        x |= (tail as u128) << (M - 1 + last_b64_bits);
        for _ in 0..(M - 1 + last_b64_bits) {
            let idx = (x as u64 & mask) as usize;
            bucket[idx] += 1;
            x >>= 1;
        }
    }

    if let Ok(mut pat) = sample.patterns.lock() {
        if pat.len() <= M {
            pat.resize(M + 1, None);
        }
        pat[M] = Some(bucket.clone());
    }

    bucket
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_popcount(e: &[u8]) -> u64 {
        e.iter()
            .map(|a| {
                let a = *a;
                let mut s = 0;
                for i in 0..8 {
                    s += (a >> i) & 1;
                }
                s as u64
            })
            .sum()
    }

    #[test]
    fn test_popcount() {
        for n in 1..10 {
            let mut e = Vec::new();
            let mut b = Vec::new();
            for i in 1010u64..(10000 + n) {
                e.push((i * i * i % 7 % 2) as u8);
                b.push((i * i * i * i) as u8);
            }
            assert_eq!(
                naive_popcount(e.as_slice()),
                popcount_epsilon(e.as_slice())
            );
            assert_eq!(naive_popcount(b.as_slice()), popcount(b.as_slice()));
        }
    }

    extern crate test;
    use test::Bencher;
    #[bench]
    fn bench_popcount(b: &mut Bencher) {
        let mut e = Vec::new();
        for i in 0u64..1000000 {
            e.push((i * i * i % 7 % 2) as u8);
        }
        // 11,890.95 ns/iter
        b.iter(|| {
            test::black_box(popcount_epsilon(e.as_slice()));
        });
    }

    #[bench]
    fn bench_popcount_dense(b: &mut Bencher) {
        let mut e = Vec::new();
        for i in 0u64..1000000 / 8 {
            e.push((i * i * i) as u8);
        }
        // 1,026.40 ns/iter ns/iter
        b.iter(|| {
            test::black_box(popcount(e.as_slice()));
        });
    }
}
