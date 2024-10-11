use std::fmt::Display;

use super::rank_u64::*;
use crate::{Sample, TestResult, USE_U8};

// The fixed 32x32 matrix dimension
const MATRIX_DIM: usize = 32;
const MATRIX_SIZE: usize = MATRIX_DIM * MATRIX_DIM;

/// 矩阵秩检测
pub(crate) fn rank(sample: &Sample) -> TestResult {
    if USE_U8 {
        rank_u8(sample)
    } else {
        rank_u64(sample)
    }
}
#[cfg(test)]
pub(crate) fn rank_epsilon(sample: &Sample) -> TestResult {
    use crate::{igamc, powi};
    let e = &sample.e;

    let N = e.len() / MATRIX_SIZE;

    if N == 0 {
        return TestResult::default();
    }

    let p_32 = 0.2888;
    let p_31 = 0.5776;
    let p_30 = 0.1336;

    let mut F_32 = 0;
    let mut F_31 = 0;
    let mut m = Matrix::default();
    for chunk in e.chunks_exact(MATRIX_SIZE) {
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

struct Matrix {
    data: [u8; MATRIX_SIZE], // the elements of the matrix
    row_index: [usize; MATRIX_DIM], // [0,32,64,...], the beginning of the row index in data
}

impl Default for Matrix {
    fn default() -> Self {
        Self {
            data: [0; MATRIX_SIZE],
            row_index: [0; MATRIX_DIM],
        }
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..MATRIX_DIM {
            write!(
                f,
                "{:?}\n",
                &self.data[self.row_index[i]..self.row_index[i] + MATRIX_DIM]
            )?;
        }
        Ok(())
    }
}

impl Matrix {
    #[inline(always)]
    fn get(&self, i: usize, j: usize) -> u8 {
        return self.data[self.row_index[i] + j];
    }

    #[inline(always)]
    fn set(&mut self, i: usize, j: usize, v: u8) {
        self.data[self.row_index[i] + j] = v;
    }

    fn from_slice(&mut self, s: &[u8]) {
        self.data.copy_from_slice(&s[..MATRIX_SIZE]);
        for i in 0..MATRIX_DIM {
            self.row_index[i] = i * MATRIX_DIM;
        }
    }

    fn swap_row(&mut self, i: usize, j: usize) {
        (self.row_index[i], self.row_index[j]) =
            (self.row_index[j], self.row_index[i]);
    }

    // add row2 to row1, starting from start_pos
    fn add_row(&mut self, row1: usize, row2: usize, start_pos: usize) {
        for k in start_pos..MATRIX_DIM {
            let v = (self.get(row1, k) + self.get(row2, k)) % 2;
            self.set(row1, k, v);
        }
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
                self.find_and_swap_row(i, j);
            }

            if self.get(i, j) == 0 {
                // move to next column
                j += 1;
                continue;
            }

            // eleminate
            for k in i + 1..MATRIX_DIM {
                if self.get(k, j) == 1 {
                    self.add_row(k, i, j); // eleminate row k
                }
            }
            i += 1;
            j += 1;
        }

        i
    }
}

#[cfg(test)]
mod test {
    use super::Matrix;
    #[test]
    fn test_matrix() {
        let mut e = Vec::new();
        for i in 0..32 {
            for _ in 0..32 {
                e.push(i as u8);
            }
        }
        let e = e;
        let mut m = Matrix::default();
        m.from_slice(e.as_slice());
        // println!("{}", m);

        m.swap_row(0, 1);
        // println!("{}", m);
        // println!("{:?}", e);
    }
}
