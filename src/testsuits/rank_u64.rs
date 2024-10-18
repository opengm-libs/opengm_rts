
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
        m.from_slice_u64(chunk);
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


/// 矩阵秩检测
pub(crate) fn rank_u8(sample: &Sample) -> TestResult {
    let b8 = &sample.b;
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
    let chunk_size = MATRIX_DIM * MATRIX_DIM / 8;
    for chunk in b8[..N*chunk_size].chunks_exact(chunk_size) {
        m.from_slice_u8(chunk);
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

#[derive(Debug, Default, Clone)]
struct Matrix {
    m: [u32; MATRIX_DIM],
}

impl From<&[u8]> for Matrix{
    fn from(value: &[u8]) -> Self {
        let mut m = [0u32;MATRIX_DIM];

        for i in 0..32{
            m[i] = u32::from_be_bytes(value[i*4..(i+1)*4].try_into().unwrap());
        }
        Matrix{m}
    }
}

impl From<&Vec<u8>> for Matrix{
    fn from(value: &Vec<u8>) -> Self {
        let mut m = [0u32;MATRIX_DIM];

        for i in 0..32{
            m[i] = u32::from_be_bytes(value[i*4..(i+1)*4].try_into().unwrap());
        }
        Matrix{m}
    }
}


impl Matrix {
    #[inline(always)]
    fn from_slice_u64(&mut self, m: &[u64]) {
        for (i, x) in m.iter().enumerate() {
            for j in 0..64 / MATRIX_DIM {
                self.m[2*i+j] = ((*x) >> (MATRIX_DIM * j)) as u32;
            }
        }
    }
    #[inline(always)]
    fn from_slice_u8(&mut self, m: &[u8]) {
        for i in 0..32{
            self.m[i] = u32::from_ne_bytes(m[i*4..(i+1)*4].try_into().unwrap());
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

    #[inline(always)]
    fn add_row(&mut self, row1: usize, row2: usize) {
        self.m[row1] ^= self.m[row2];
    }

    // assume [i,j] = 0
    // find a nonzero row k > i which at column index j and swap with j.
    // if find and swap such row, return true.
    #[inline(always)]
    fn find_and_swap_row(&mut self, i: usize, maskj: u32) -> bool {
        for k in (i + 1)..MATRIX_DIM {
            if self.m[k] & maskj != 0{
                self.swap_row(k, i);
                return true;
            }
        }
        false
    }
    #[inline(always)]
    fn compute_rank(&mut self) -> usize {
        let mut i = 0;
        let mut j = 0;
        while j < MATRIX_DIM {
            let mask = 1<<(31-j);
            if self.m[i] & mask == 0{
                if !self.find_and_swap_row(i, mask){
                    // move to next column
                    j += 1;
                    continue;
                }
            }

            // eleminate
            for k in i + 1..MATRIX_DIM {
                if self.m[k] & mask != 0{
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
        assert_eq!(tv1.1, rank_u8(&sample));
    }


    #[test]
    fn test_equal() {
        for nbits in 250..1000 {
            let sample: Sample = E[..nbits * 8].into();
            assert_eq!(rank_epsilon(&sample), rank_u64(&sample));
            assert_eq!(rank_u8(&sample), rank_u64(&sample));
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
    fn bench_rank_u64(b: &mut Bencher) {
        let sample: Sample = E.into();
        // let m = sample.b.into();

        // 6,719.88 ns/iter
        b.iter(|| {
            test::black_box(rank_u64(&sample));
        });
    }

    #[bench]
    fn bench_rank_u8(b: &mut Bencher) {
        let sample: Sample = E.into();

        // 1,889,683.30
        b.iter(|| {
            test::black_box(rank_u8(&sample));
        });
    }

    #[bench]
    fn bench_rank(b: &mut Bencher) {
        let sample: Sample = E.into();
        let m0:Matrix = (&sample.b).into();
        let mut m = Matrix::default();

        // 392.75 ns/iter
        b.iter(|| {
            m.m.copy_from_slice(&m0.m);
            test::black_box(m.compute_rank());
        });
    }

}