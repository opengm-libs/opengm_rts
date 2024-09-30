
use super::util::*;
use crate::{Sample, TestResult};

// The fixed 32x32 matrix dimension
const MATRIX_DIM: usize = 32;
const MATRIX_SIZE: usize = MATRIX_DIM * MATRIX_DIM;

/// 矩阵秩检测
pub(crate) fn rank_u64(sample: &Sample) -> TestResult {
    let b64 = &sample.b64;
    let n = sample.len();

    let N = n / (MATRIX_DIM * MATRIX_DIM);
    
    if N == 0 {
        return TestResult::default();
    }
    
    let p_32 = 0.2888;
    let p_31 = 0.5776;
    let p_30 = 0.1336;
    
    let mut F_32 = 0;
    let mut F_31 = 0;
    let mut m = Matrix::default();
    let chunk_size = MATRIX_DIM * MATRIX_DIM / 64;
    for chunk in b64[..N*chunk_size].chunks_exact(chunk_size) {
        m.from_slice(chunk);
        let rank = m.compute_rank();
        if 32 == rank {
            F_32 += 1
        }
        if 31 == rank {
            F_31 += 1
        }
    }
    let F_30 = N - F_32 - F_31;

    let F_32 = F_32 as f64;
    let F_31 = F_31 as f64;
    let F_30 = F_30 as f64;
    let N = N as f64;

    let mut v = 0.0;
    v += powi(F_32 - N * p_32, 2) / (N * p_32);
    v += powi(F_31 - N * p_31, 2) / (N * p_31);
    v += powi(F_30 - N * p_30, 2) / (N * p_30);

    let pv = igamc(1.0, v / 2.0);
    TestResult { pv, qv: pv }
}

#[derive(Debug, Default)]
struct Matrix {
    m: [u32; MATRIX_DIM],
}

impl Matrix {
    #[inline(always)]
    fn from_slice(&mut self, m: &[u64]) {
        for (i, x) in m.iter().enumerate() {
            for j in 0..64 / MATRIX_DIM {
                self.m[2*i+j] = ((*x) >> (MATRIX_DIM * j)) as u32;
            }
        }
    }

    #[inline(always)]
    fn get(&self, i: usize, j: usize) -> u8 {
        ((self.m[i] >> (31-j)) & 1) as u8
    }

    // v = 0,1.
    #[inline(always)]
    fn set(&mut self, i: usize, j: usize, v: u8) {
        self.m[i] &= (!1)<<j; // clear the (i,j) bit
        self.m[i] |= (v as u32) << (31-j);
    }

    #[inline(always)]
    fn swap_row(&mut self, i: usize, j: usize) {
        (self.m[i], self.m[j]) = (self.m[j], self.m[i])
    }

    // add row2 to row1, starting from start_pos
    // #[inline(always)]
    // fn add_row(&mut self, row1: usize, row2: usize, start_pos: usize) {
    //     // for k in start_pos..MATRIX_DIM {
    //     //     let v = (self.get(row1, k) + self.get(row2, k)) % 2;
    //     //     self.set(row1, k, v);
    //     // }

    //     // 0..01111 with start_pos '0'.
    //     let mask = (!0)>>start_pos;
    //     let mut row = self.m[row1] ^ self.m[row2];
    //     self.m[row1] &= !mask;
    //     self.m[row1] |= row & mask;
    // }

    #[inline(always)]
    fn add_row(&mut self, row1: usize, row2: usize) {
        self.m[row1] ^= self.m[row2];
    }

    // assume [i,j] = 0
    // find a nonzero row k > i which at column index j and swap with j.
    // if find and swap such row, return true.
    fn find_and_swap_row(&mut self, i: usize, j: usize) -> bool {
        for k in (i + 1)..MATRIX_DIM {
            if self.get(k, j) != 0 {
                self.swap_row(k, i);
                return true;
            }
        }
        false
    }

    fn compute_rank(&mut self) -> usize {
        let mut i = 0;
        let mut j = 0;
        while j < MATRIX_DIM {
            if self.get(i, j) == 0 {
                if !self.find_and_swap_row(i, j){
                    // move to next column
                    j += 1;
                    continue;
                }
            }

            // eleminate
            for k in i + 1..MATRIX_DIM {
                if self.get(k, j) == 1 {
                    self.add_row(k, i); // eleminate row k
                }
            }
            i += 1;
            j += 1;
        }

        i
    }
}

#[cfg(test)]
mod tests {
    use super::super::rank::*;
    use super::{*, super::tests::*};
    use crate::{test_data::E, Sample};

    #[test]
    fn test_basic() {
        let tv = get_test_vec(crate::TestFuncs::Rank);
        let sample: Sample = tv.0.into();
        assert_eq!(rank_epsilon(&sample), tv.2);
        assert_eq!(rank_u64(&sample), tv.2);
    }

    #[test]
    fn test_e() {
        let tv1 = get_test_vec_e(crate::TestFuncs::Rank);
        let sample: Sample = E.into();
        assert_eq!(tv1.1, rank_epsilon(&sample));
        assert_eq!(tv1.1, rank_u64(&sample));
    }

    #[test]
    fn test_equal() {
        for nbits in 250..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(rank_epsilon(&sample), rank_u64(&sample));
        }
    }
}
