use super::util::*;
use crate::{Sample, TestResult};

/// 近似熵检测
pub(crate) fn approximate_entropy(sample: &Sample, m:i32) -> TestResult {
    let n = sample.e.len();
	let apen = phi(sample.e.as_slice(), m)-phi(sample.e.as_slice(), m+1);
	
	let v = 2.0 * (n as f64) * (ln(2.0) - apen);
    let pv = igamc(powi(2.0, m-1), v / 2.0);
    TestResult {
        pv,
        qv: pv,
    }
}

// m = 2, 5
const MAX_M: usize = 6;
fn phi(e: &[u8], m: i32) -> f64{
    if m == 0  {
        return 0.0;
    }
    if m as usize > MAX_M {
        panic!("Invalid m")
    }

    let n = e.len();
    let m = m as usize;
    let mut probs = [0usize; 1<<MAX_M];

    // buf = {e[n-m+1], ..., e[n-1], e[0], e[1], e[m-2]}
    let mut buf = [0u8; 2*MAX_M];    // buf.len = m-1 + m-1, only 2m-2 are need.
    for i in 0..(2*m-2){
        buf[i] = e[(n-m+1+i) % n];
    }

    for w in e.windows(m){
        probs[get_partten(w)] += 1;
    }

    for w in buf[..(2*m-2)].windows(m){
        probs[get_partten(w)] += 1;
    }

    let mut sum = 0.0;
    for i in 0..(1<<m) {
        if probs[i] == 0{
            continue;
        }
        let c = probs[i] as f64 / n as f64;
        sum += c * ln(c);
    }
    sum
}

fn get_partten(e: &[u8]) -> usize {
    let mut k = 0;
    for i in e {
        k <<= 1;
        k |= *i as usize;
    }
    k
}