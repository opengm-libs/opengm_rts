use super::util::*;
use crate::{Sample, TestResult};

// Note: for k = 2^n - 1, ei = epsilon[i]^epsilon[i+1]^epsilon[i+2]^... ^epsilon[i+k]
//  0: 0        1        2         3          4          5         6  ...
//  1: 01       12       23        34         45         56        67 ...
//  2: 02       13       24        35         46         57        68 ... 
//  3: 0123     1234     2345      3456       ...
//  4: 04       15       26        37         ...
//  5: 0145     1256     2367      3478       ...
//  6: 0246     1357     2468      3579       ...
//  7: 01234567 12345678 ...
//  8: 08       19       2(10)     3(11)      ...       
//  9: 0189     129(10)
// 10: 028(10)  ...
// 11: 
// 12: 
// 13: 
// 14: 
// 15: 0-15     1-16

/// 二元推导检测
pub(crate) fn binary_derivative(sample: &Sample, k:i32)-> TestResult {
    let epsilon = &sample.e;
    let n = epsilon.len();
    let k = k as usize;
    let pv;
    let qv;

    // k is of the form 2^m - 1.
    if 1 == (k+1).count_ones(){
        // e'[i] = epsilon[i]^epsilon[i+1]^epsilon[i+2]^... ^epsilon[i+7]
        let mut ei = 0u8;
        // epsilon[0] ^ epsilon[1] ^ ... ^ epsilon[k];
        for i in 0..(k+1) {
            ei ^= epsilon[i];
        }
        let mut sum = ei as i64;
        for i in 1..(n - k){
            ei = ei ^ epsilon[i - 1] ^ epsilon[i + k];
            sum += ei as i64;
        }
        sum = 2 * sum - (n - k) as i64;
        pv = erfc(abs(sum as f64) / sqrt((n - k) as f64) / SQRT2);
        qv = erfc((sum as f64) / sqrt((n - k) as f64) / SQRT2)/2.0;
    }else{
        // make a copy
        let mut tmp_epsilon = epsilon.clone();

        for j in 1..(k+1){
            for i in 0..(n-j){
                tmp_epsilon[i] ^= tmp_epsilon[i + 1];
            }
        }
        let mut sum = 0.0;
        for i in 0..(n-k){
            sum += (2 * (tmp_epsilon[i] as i8) - 1) as f64;
        }

        pv = erfc(abs(sum) / sqrt((n - k) as f64) / SQRT2);
        qv = erfc(sum / sqrt((n - k) as f64) / SQRT2)/2.0;
    }

    TestResult {
        pv1: pv,
        qv1: qv,
        pv2: None,
        qv2: None,
    }

}