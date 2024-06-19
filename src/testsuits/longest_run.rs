use super::util::*;
use crate::Sample;

#[inline(always)]
fn longest_run_internal(e: &[u8], K: i32, M: i32, V: &[i32], pi: &[f64]) -> (f64, f64) {
    let n = e.len();
    let N = n / (M as usize);
    let mut nu0 = [0; 7]; // K at most 6
    let mut nu1 = [0; 7]; // K at most 6

    for i in 0..N {
        let block = &e[i * (M as usize)..];
        let mut longestRun0 = 0;
        let mut longestRun1 = 0;
        let mut currentRun0 = 0;
        let mut currentRun1 = 0;
        for j in 0..(M as usize) {
            if block[j] == 0 {
                currentRun0 += 1;
                longestRun0 = saturating(longestRun0, currentRun0, longestRun0);
                currentRun1 = 0
            } else {
                currentRun1 += 1;
                longestRun1 = saturating(longestRun1, currentRun1, longestRun1);
                currentRun0 = 0;
            }
        }
        nu0[(saturating(longestRun0, V[0], V[K as usize]) - V[0]) as usize] += 1;
        nu1[(saturating(longestRun1, V[0], V[K as usize]) - V[0]) as usize] += 1;
    }

    let mut chi2 = [0.0,0.0];
    for i in 0..(K+1) {
        let N = N as f64;
        let nu0i = nu0[i as usize] as f64;
        let nu1i = nu1[i as usize] as f64;
        let pii = pi[i as usize];
        chi2[0] += (nu0i - N * pii) * (nu0i - N * pii) / (N * pii);
        chi2[1] += (nu1i - N * pii) * (nu1i - N * pii) / (N * pii);
    }

    (igamc(K as f64 / 2.0, chi2[0] / 2.0),igamc(K as f64 / 2.0, chi2[1] / 2.0))

}

// 7. 块内最大游程检测
pub(crate) fn longest_run(sample: &Sample) -> (f64, f64) {
    if sample.e.len() < 128 {
        panic!("longest run test: n too short, at least 128\n");
    }

    let n = sample.e.len();
    let K;
    let M;
    let mut pi = [0.0; 7];
    let mut V = [0; 7];

    // ref: B.7 of GM/T 0005-2021
    if n < 6272 {
        K = 3;
        M = 8;
        V[0] = 1;
        V[1] = 2;
        V[2] = 3;
        V[3] = 4;
        pi[0] = 0.2148;
        pi[1] = 0.3672;
        pi[2] = 0.2305;
        pi[3] = 0.1875;
    } else if n < 750000 {
        K = 5;
        M = 128;
        V[0] = 4;
        V[1] = 5;
        V[2] = 6;
        V[3] = 7;
        V[4] = 8;
        V[5] = 9;
        pi[0] = 0.1174;
        pi[1] = 0.2430;
        pi[2] = 0.2494;
        pi[3] = 0.1752;
        pi[4] = 0.1027;
        pi[5] = 0.1124;
    } else {
        K = 6;
        M = 10000;
        V[0] = 10;
        V[1] = 11;
        V[2] = 12;
        V[3] = 13;
        V[4] = 14;
        V[5] = 15;
        V[6] = 16;
        pi[0] = 0.086632;
        pi[1] = 0.208201;
        pi[2] = 0.248419;
        pi[3] = 0.193913;
        pi[4] = 0.121458;
        pi[5] = 0.068011;
        pi[6] = 0.073366;
    }

    longest_run_internal(&sample.e, K, M, &V, &pi)
}
