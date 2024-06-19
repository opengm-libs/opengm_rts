mod math;

pub(crate) use math::*;

// returns the sum of epsilon[0] + ... + epsilon[n-1]
// assumes that epsilon[i] = 0 or 1.
pub(crate) fn popcount(e: &[u8]) -> u64 {
    let mut sum = 0u64;

    let (bigchunks, remainder) = e.as_chunks::<2040>();
    for bigchunk in bigchunks {
        let mut tmp_sum = 0u64;
        for chunk in unsafe { bigchunk.as_chunks_unchecked::<8>() } {
            tmp_sum += u64::from_ne_bytes(*chunk)
        }
        for i in 0..8 {
            sum += (tmp_sum >> 8 * i) & 0xff;
        }
    }

    let (chunks, remainder) = remainder.as_chunks::<8>();
    let mut tmp_sum = 0u64;
    for chunk in chunks {
        tmp_sum += u64::from_ne_bytes(*chunk)
    }
    for i in 0..8 {
        sum += (tmp_sum >> 8 * i) & 0xff;
    }

    for a in remainder {
        sum += *a as u64;
    }

    sum
}



#[inline(always)]
pub(crate) fn saturating<T: core::cmp::PartialOrd>(n: T, min: T, max: T) -> T {
    if n > max {
        max
    } else if n < min {
        min
    } else {
        n
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_popcount(e: &[u8]) -> u64 {
        e.iter().map(|a| *a as u64).sum()
    }

    #[test]
    fn test_popcount() {
        for n in 1..10 {
            let mut e = Vec::new();
            for i in 0u64..(10000 + n) {
                e.push((i * i * i % 7 % 2) as u8)
            }
            assert_eq!(naive_popcount(e.as_slice()), popcount(e.as_slice()));
        }
    }

    extern crate test;
    use test::Bencher;
    #[bench]
    fn bench_popcount(b: &mut Bencher) {
        let mut e = Vec::new();
        for i in 0u64..(1000000) {
            e.push((i * i * i % 7 % 2) as u8)
        }
        b.iter(|| popcount(e.as_slice()));
    }

    #[bench]
    fn bench_naive_popcount(b: &mut Bencher) {
        let mut e = Vec::new();
        for i in 0u64..(1000000) {
            e.push((i * i * i % 7 % 2) as u8)
        }
        b.iter(|| naive_popcount(e.as_slice()));
    }
}
