mod math;

pub(crate) use math::*;

// returns the sum of epsilon[0] + ... + epsilon[n-1]
// assumes that epsilon[i] = 0 or 1.
pub(crate) fn popcount(e: &[u8]) -> u64 {
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
}
