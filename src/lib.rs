#![allow(dead_code)]
#![allow(non_snake_case)]
#![feature(test)]
#![feature(slice_as_chunks)]
#![feature(float_gamma)]
#![feature(new_uninit)]

extern crate core;
use testsuits::util::popcount;
mod testsuits;

enum TestFuncs {
    Frequency,          // 1.frequency 频度检测
    BlockFrequency,     // 2.block_frequency 块内频度检测
    Poker,              // 3.poker 扑克检测
    Serial,             // 4.serial
    Runs,               // 5.runs 游程分布检测
    RunsDistribution,   // 6.runs_distribution 游程分布检测
    LongestRun,         // 7.longest_run 块内最大游程检测
    BinaryDerivative,   // 8.binary_derivative 二元推导检测
    Autocorrelation,    // 9.autocorrelation 自相关检测
    Rank,               // 10.rank 矩阵秩检测
    CumulativeSums,     // 11.cumulative_sums 累加和检测
    ApproximateEntropy, // 12.approximate_entropy 近似熵检测
    LinearComplexity,   // 13.linear_complexity 线性复杂度检测
    Universal,          // 14.universal Maurer通用统计检测
    DiscreteFourier,    // 15.discrete_fourier 离散傅立叶检测
}

pub struct Sample {
    e: Vec<u8>,
    pop: u64, // popcount for e
}

// from bit string("10100...")
impl From<&str> for Sample {
    fn from(bit_string: &str) -> Self {
        let mut e = Vec::with_capacity(bit_string.len());
        for c in bit_string.chars() {
            e.push(c as u8 - '0' as u8);
        }
        let pop = popcount(&e.as_slice());
        Sample { e: e, pop: pop }
    }
}

impl From<&[u8]> for Sample {
    fn from(byte_slice: &[u8]) -> Self {
        let mut e = Vec::with_capacity(byte_slice.len() * 8);
        for b in byte_slice {
            for j in 0..8 {
                e.push((b >> (7 - j)) & 1)
            }
        }
        let pop = popcount(&e.as_slice());
        Sample { e: e, pop: pop }
    }
}

impl Sample {
    pub fn new(bytes: &[u8]) -> Sample {
        return bytes.into();
    }

    pub fn frequency(&self) -> f64 {
        testsuits::frequency(self)
    }
    pub fn block_frequency(&self, m: i32) -> f64 {
        testsuits::block_frequency(self, m)
    }

    pub fn poker(&self, m: i32) -> f64 {
        testsuits::poker(self, m)
    }

    pub fn serial(&self, m: i32) -> (f64, f64) {
        testsuits::serial(self, m)
    }

    pub fn runs(&self) -> f64 {
        testsuits::runs(self)
    }

    pub fn runs_distribution(&self) -> f64 {
        testsuits::runs_distribution(self)
    }

    pub fn longest_run(&self) -> (f64, f64) {
        testsuits::longest_run(self)
    }

    pub fn binary_derivative(&self, k: i32) -> f64 {
        testsuits::binary_derivative(self, k)
    }

    pub fn autocorrelation(&self, d: i32) -> f64 {
        testsuits::autocorrelation(self, d)
    }

    pub fn rank(&self) -> f64 {
        testsuits::rank(self)
    }

    pub fn cumulative_sums(&self) -> (f64, f64) {
        testsuits::cumulative_sums(self)
    }

    pub fn approximate_entropy(&self, m: i32) -> f64 {
        testsuits::approximate_entropy(self, m)
    }

    pub fn linear_complexity(&self, m: i32) -> f64 {
        testsuits::linear_complexity(self, m)
    }

    pub fn universal(&self) -> f64 {
        testsuits::universal(self)
    }

    pub fn discrete_fourier(&self) -> f64 {
        testsuits::discrete_fourier(self)
    }
}

#[cfg(test)]
mod test_data;

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::test_data::E;
    use super::*;

    struct TestVec {
        epsilon_str: &'static str,
        param: i32,
        pvalue: [f64; 2],
    }

    const TEST_VEC:&[TestVec] = &[
        // 1.frequency 频度检测
        TestVec{
            epsilon_str:"11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param:0,
            pvalue:[0.215925,0.0]},
        // 2.block_frequency 块内频度检测
        TestVec{
            epsilon_str : "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            param : 10,
            pvalue : [0.706438,0.0]},
        // 3.poker 扑克检测
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param : 4,
            pvalue : [0.213734,0.0]},
        // 4.serial
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param : 2,
            pvalue : [0.436868, 0.723674]},
        // 5.runs 游程总数检测
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param:0,
            pvalue : [0.620729,0.0]},
        // 6.runs_distribution 游程分布检测
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param:0,
            pvalue : [0.970152,0.0]},
        // 7.longest_run 块内最大游程检测
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param:0,
            pvalue : [0.839299,0.180598]},// pvalue[0] = 块内最大“0”游程检测, pvalue[1] = 块内最大“1”游程检测
        // 8.binary_derivative 二元推导检测
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param : 3,
            pvalue : [0.039669,0.0]},
        // 9.auto_correlation 自相关检测
        TestVec{
            epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            param : 1,
            pvalue : [0.790080,0.0]},
        // 10.rank 矩阵秩检测
        TestVec{
            epsilon_str : E,
            param: 1000000,
            pvalue : [0.307543,0.0]},
        // 11.cumulative_sums 累加和检测
        TestVec{
            epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            param:0,
            pvalue : [0.219194, 0.114866]},// pvalue[0] = 前向累加和检测, pvalue[1] = 后向累加和检测
        // 12.approximate_entropy 近似熵检测
        TestVec{
            epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            param : 2,
            pvalue : [0.235301,0.0],},
        // 13.linear_complexity 线性复杂度检测
        TestVec{
            epsilon_str : E,
            param : 1000,
            pvalue : [0.844721,0.0]},
        // 14.universal Maurer通用统计检测
        TestVec{
            epsilon_str : E,
            param:0,
            pvalue : [0.282568,0.0]},
        // 15.discrete_fourier_transform 离散傅立叶检测
        TestVec{
            epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            param:0,
            pvalue : [0.654721,0.0]}];

    fn assert_eq_f64(a: f64, b: f64) {
        if (a - b).abs() > 0.001 {
            panic!("{} != {}", a, b);
        }
    }

    #[test]
    fn test_frequency() {
        let tv = &TEST_VEC[TestFuncs::Frequency as usize];
        let sample: Sample = tv.epsilon_str.into();
        assert_eq_f64(sample.frequency(), tv.pvalue[0]);
    }

    #[test]
    fn test_block_frequency() {
        let tv = &TEST_VEC[TestFuncs::BlockFrequency as usize];
        let sample: Sample = tv.epsilon_str.into();
        assert_eq_f64(sample.block_frequency(tv.param), tv.pvalue[0]);
    }

    #[test]
    fn test_poker() {
        let tv = &TEST_VEC[TestFuncs::Poker as usize];
        let sample: Sample = tv.epsilon_str.into();
        assert_eq_f64(sample.poker(tv.param), tv.pvalue[0]);
    }

    #[test]
    fn test_serial() {
        let tv = &TEST_VEC[TestFuncs::Serial as usize];
        let sample: Sample = tv.epsilon_str.into();
        let (pvalue1, pvalue2) = sample.serial(tv.param);
        assert_eq_f64(pvalue1, tv.pvalue[0]);
        assert_eq_f64(pvalue2, tv.pvalue[1]);
    }

    #[test]
    fn test_runs() {
        let tv = &TEST_VEC[TestFuncs::Runs as usize];
        let sample: Sample = tv.epsilon_str.into();
        assert_eq_f64(sample.runs(), tv.pvalue[0]);
    }

    #[test]
    fn test_runs_distribution() {
        let tv = &TEST_VEC[TestFuncs::RunsDistribution as usize];
        let sample: Sample = tv.epsilon_str.into();
        assert_eq_f64(sample.runs_distribution(), tv.pvalue[0]);
    }

    #[test]
    fn test_longest_run() {
        let tv = &TEST_VEC[TestFuncs::LongestRun as usize];
        let sample: Sample = tv.epsilon_str.into();
        let (pv0, pv1) = sample.longest_run();
        assert_eq_f64(pv0, tv.pvalue[0]);
        assert_eq_f64(pv1, tv.pvalue[1]);
    }

    #[test]
    fn test_binary_derivative() {
        let tv = &TEST_VEC[TestFuncs::BinaryDerivative as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.binary_derivative(tv.param);
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    #[test]
    fn test_autocorrelation() {
        let tv = &TEST_VEC[TestFuncs::Autocorrelation as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.autocorrelation(tv.param);
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    #[test]
    fn test_rank() {
        let tv = &TEST_VEC[TestFuncs::Rank as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.rank();
        // let pv = sample.c_rank();
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    #[test]
    fn test_cumulative_sums() {
        let tv = &TEST_VEC[TestFuncs::CumulativeSums as usize];
        let sample: Sample = tv.epsilon_str.into();
        let (pv0, pv1) = sample.cumulative_sums();
        assert_eq_f64(pv0, tv.pvalue[0]);
        assert_eq_f64(pv1, tv.pvalue[1]);
    }

    #[test]
    fn test_approximate_entropy() {
        let tv = &TEST_VEC[TestFuncs::ApproximateEntropy as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.approximate_entropy(tv.param);
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    #[test]
    fn test_linear_complexity() {
        let tv = &TEST_VEC[TestFuncs::LinearComplexity as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.linear_complexity(tv.param);
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    #[test]
    fn test_universal() {
        let tv = &TEST_VEC[TestFuncs::Universal as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.universal();
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    #[test]
    fn test_discrete_fourier() {
        let tv = &TEST_VEC[TestFuncs::DiscreteFourier as usize];
        let sample: Sample = tv.epsilon_str.into();
        let pv = sample.discrete_fourier();
        assert_eq_f64(pv, tv.pvalue[0]);
    }

    extern crate test;
    use test::Bencher;
    #[bench]
    fn bench_linear_complexity(b: &mut Bencher) {
        let sample: Sample = E.into();
        b.iter(|| sample.linear_complexity(1000));
    }

    #[bench]
    fn bench_rank(b: &mut Bencher) {
        let sample: Sample = E.into();
        b.iter(|| sample.rank());
    }

    #[bench]
    fn bench_discrete_fourier(b: &mut Bencher) {
        let sample: Sample = E.into();
        b.iter(|| sample.discrete_fourier());
    }

    // Test for all 15 test functions.
    #[test]
    fn test_all() {
        let sample: Sample = E.into();
        let mut pv = Vec::new();

        let start = Instant::now();
        let last = start;

        pv.clear();
        pv.push(sample.frequency());
        let now = Instant::now();
        println!("frequency: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.block_frequency(10000));
        let now = Instant::now();
        println!("block_frequency: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.poker(8));
        let now = Instant::now();
        println!("poker: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.serial(5).0);
        let now = Instant::now();
        println!("serial: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.runs());
        let now = Instant::now();
        println!("runs: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.runs_distribution());
        let now = Instant::now();
        println!("runs_distribution: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.longest_run().0);
        let now = Instant::now();
        println!("longest_run: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.binary_derivative(7));
        let now = Instant::now();
        println!("binary_derivative: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.autocorrelation(16));
        let now = Instant::now();
        println!("autocorrelation: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.rank());
        let now = Instant::now();
        println!("rank: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.cumulative_sums().0);
        let now = Instant::now();
        println!("cumulative_sums: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.approximate_entropy(5));
        let now = Instant::now();
        println!("approximate_entropy: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.linear_complexity(1000));
        let now = Instant::now();
        println!("linear_complexity: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.universal());
        let now = Instant::now();
        println!("universal: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.discrete_fourier());
        let now = Instant::now();
        println!("discrete_fourier: {:.2} s", (now - last).as_secs_f64());
        
        println!("total: {:.2} s", (Instant::now() - start).as_secs_f64());
    }
}
