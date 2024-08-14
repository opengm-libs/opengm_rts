#![allow(dead_code)]
#![allow(non_snake_case)]

mod tester;
mod testsuits;
pub use tester::*;
use testsuits::util::popcount;
use testsuits::util::*;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum TestFuncs {
    Frequency,              // 1. 频度检测
    BlockFrequency,         // 2. 块内频度检测
    Poker,                  // 3. 扑克检测
    Serial1,                // 4. 重叠子序列p_value1检测
    Serial2,                // 4. 重叠子序列p_value2检测
    Runs,                   // 5. 游程分布检测
    RunsDistribution,       // 6. 游程分布检测
    LongestRun0,            // 7. 块内最大0游程检测
    LongestRun1,            // 7. 块内最大1游程检测
    BinaryDerivative,       // 8. 二元推导检测
    Autocorrelation,        // 9. 自相关检测
    Rank,                   // 10. 矩阵秩检测
    CumulativeSumsForward,  // 11. 前向累加和检测
    CumulativeSumsBackward, // 11. 后向累加和检测
    ApproximateEntropy,     // 12. 近似熵检测
    LinearComplexity,       // 13. 线性复杂度检测
    Universal,              // 14. Maurer通用统计检测
    DiscreteFourier,        // 15. 离散傅立叶检测
}

#[derive(Default, Debug)]
pub struct TestResult {
    pub pv: f64,
    pub qv: f64,
}

impl TestResult {
    pub fn pass(&self, alpha: f64) -> bool {
        self.pv >= alpha
    }
}

pub struct Sample {
    e: Vec<u8>,
    
    // popcount for e
    pop: u64,
}

// from bit string("10100...")
impl From<&str> for Sample {
    fn from(bit_string: &str) -> Self {
        let mut e = Vec::with_capacity(bit_string.len());
        for c in bit_string.chars() {
            e.push(c as u8 - '0' as u8);
        }
        let pop = popcount(&e.as_slice());
        Sample {
            e: e,
            pop: pop,
        }
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
        Sample {
            e: e,
            pop: pop,
        }
    }
}

impl Sample {
    pub fn new(bytes: &[u8]) -> Sample {
        return bytes.into();
    }

    pub fn bits(&self) -> usize {
        return self.e.len();
    }
    pub fn frequency(&self) -> TestResult {
        testsuits::frequency(self)
    }

    pub fn block_frequency(&self, m: i32) -> TestResult {
        testsuits::block_frequency(self, m)
    }

    pub fn poker(&self, m: i32) -> TestResult {
        testsuits::poker(self, m)
    }

    pub fn serial1(&self, m: i32) -> TestResult {
        testsuits::serial1(self, m)
    }
    pub fn serial2(&self, m: i32) -> TestResult {
        testsuits::serial2(self, m)
    }

    pub fn runs(&self) -> TestResult {
        testsuits::runs(self)
    }

    pub fn runs_distribution(&self) -> TestResult {
        testsuits::runs_distribution(self)
    }

    pub fn longest_run0(&self) -> TestResult {
        testsuits::longest_run0(self)
    }

    pub fn longest_run1(&self) -> TestResult {
        testsuits::longest_run1(self)
    }

    pub fn binary_derivative(&self, k: i32) -> TestResult {
        testsuits::binary_derivative(self, k)
    }

    pub fn autocorrelation(&self, d: i32) -> TestResult {
        testsuits::autocorrelation(self, d)
    }

    pub fn rank(&self) -> TestResult {
        testsuits::rank(self)
    }

    pub fn cumulative_sums_forward(&self) -> TestResult {
        testsuits::cumulative_sums_forward(self)
    }
    pub fn cumulative_sums_backward(&self) -> TestResult {
        testsuits::cumulative_sums_backward(self)
    }

    pub fn approximate_entropy(&self, m: i32) -> TestResult {
        testsuits::approximate_entropy(self, m)
    }

    pub fn linear_complexity(&self, m: i32) -> TestResult {
        testsuits::linear_complexity(self, m)
    }

    pub fn universal(&self) -> TestResult {
        testsuits::universal(self)
    }

    pub fn discrete_fourier(&self) -> TestResult {
        testsuits::discrete_fourier(self)
    }
}

/// GM/T 0005-2021, 6.3 样本分布均匀性判定
/// The value of qv[i] should distribute formally in [0,1].
pub fn qvalue_distribution(qv: &[f64], k: usize) -> f64 {
    if qv.len() <= 1 {
        return 1.0;
    }

    let s = qv.len();
    let mut F = vec![0.0; k];

    // i/k <= v < (i+1)/k => i <= kv < i+1
    for v in qv {
        F[(*v * k as f64).floor() as usize] += 1.0;
    }
    let sk = s as f64 / k as f64;
    let V: f64 = F.iter().map(|f| (*f - sk) * (*f - sk) / sk).sum();
    igamc((k - 1) as f64 / 2.0, V / 2.0)
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
        pvalue: f64,
        qvalue: f64,
    }

    fn get_test_vec(f: TestFuncs) -> TestVec {
        match f{
            // 1.frequency 频度检测
            TestFuncs::Frequency => TestVec{
                epsilon_str:"11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param:0,
                pvalue: 0.215925,
                qvalue: 0.892038,
            },
            // 2.block_frequency 块内频度检测
            TestFuncs::BlockFrequency => TestVec{
                epsilon_str : "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
                param : 10,
                pvalue : 0.706438,
                qvalue: 0.706438,
            },
            // 3.poker 扑克检测
            TestFuncs::Poker => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param : 4,
                pvalue : 0.213734,
                qvalue: 0.213734,
            },
            // 4.serial
            TestFuncs::Serial1 => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param : 2,
                pvalue : 0.436868,
                qvalue:  0.436868,
            },
            // 4.serial
            TestFuncs::Serial2 => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param : 2,
                pvalue : 0.723674,
                qvalue: 0.723674,
            },
            // 5.runs 游程总数检测
            TestFuncs::Runs => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param:0,
                pvalue : 0.620729,
                qvalue: 0.310364},
            // 6.runs_distribution 游程分布检测
            TestFuncs::RunsDistribution => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param:0,
                pvalue : 0.970152,
                qvalue: 0.970152,
            },
            // 7.longest_run 块内最大0游程检测
            TestFuncs::LongestRun0 => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param:0,
                pvalue : 0.839299,
                qvalue: 0.839299},
            // 7.longest_run 块内最大1游程检测
            TestFuncs::LongestRun1 => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param:0,
                pvalue : 0.180598,
                qvalue: 0.180598,
            },
            // 8.binary_derivative 二元推导检测
            TestFuncs::BinaryDerivative => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param : 3,
                pvalue : 0.039669,
                qvalue: 0.980166},
            // 9.auto_correlation 自相关检测
            TestFuncs::Autocorrelation => TestVec{
                epsilon_str : "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
                param : 1,
                pvalue : 0.790080,
                qvalue: 0.395040,
            },
            // 10.rank 矩阵秩检测
            TestFuncs::Rank => TestVec{
                epsilon_str : E,
                param: 1000000,
                pvalue : 0.307543,
                qvalue: 0.307543,
            },
            // 11.cumulative_sums 前向累加和检测
            TestFuncs::CumulativeSumsForward => TestVec{
                epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
                param:0,
                pvalue : 0.219194,
                qvalue:  0.219194,
            },
            // 11.cumulative_sums 后向累加和检测
            TestFuncs::CumulativeSumsBackward => TestVec{
                epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
                param:0,
                pvalue : 0.114866,
                qvalue:  0.114866,
            },
            // 12.approximate_entropy 近似熵检测
            TestFuncs::ApproximateEntropy => TestVec{
                epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
                param : 2,
                pvalue : 0.235301,
                qvalue: 0.235301,
            },
            // 13.linear_complexity 线性复杂度检测
            TestFuncs::LinearComplexity => TestVec{
                epsilon_str : E,
                param : 1000,
                pvalue : 0.844721,
                qvalue: 0.844721,
            },
            // 14.universal Maurer通用统计检测
            TestFuncs::Universal => TestVec{
                epsilon_str : E,
                param:0,
                pvalue : 0.282568,
                qvalue: 0.141284,
            },
            // 15.discrete_fourier_transform 离散傅立叶检测
            TestFuncs::DiscreteFourier => TestVec{
                epsilon_str: "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
                param:0,
                pvalue: 0.654721,
                qvalue: 0.327360,
            },
        }
    }

    fn assert_eq_f64(a: f64, b: f64) {
        if (a - b).abs() > 0.001 {
            panic!("{} != {}", a, b);
        }
    }

    #[test]
    fn test_frequency() {
        let tv = get_test_vec(TestFuncs::Frequency);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.frequency();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_block_frequency() {
        let tv = get_test_vec(TestFuncs::BlockFrequency);
        let sample: Sample = tv.epsilon_str.into();
        assert_eq_f64(sample.block_frequency(tv.param).pv, tv.pvalue);
        assert_eq_f64(sample.block_frequency(tv.param).qv, tv.qvalue);
    }

    #[test]
    fn test_poker() {
        let tv = get_test_vec(TestFuncs::Poker);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.poker(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_serial1() {
        let tv = get_test_vec(TestFuncs::Serial1);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.serial1(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_serial2() {
        let tv = get_test_vec(TestFuncs::Serial2);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.serial2(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_runs() {
        let tv = get_test_vec(TestFuncs::Runs);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.runs();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_runs_distribution() {
        let tv = get_test_vec(TestFuncs::RunsDistribution);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.runs_distribution();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_longest_run0() {
        let tv = get_test_vec(TestFuncs::LongestRun0);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.longest_run0();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }
    #[test]
    fn test_longest_run1() {
        let tv = get_test_vec(TestFuncs::LongestRun1);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.longest_run1();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_binary_derivative() {
        let tv = get_test_vec(TestFuncs::BinaryDerivative);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.binary_derivative(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_autocorrelation() {
        let tv = get_test_vec(TestFuncs::Autocorrelation);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.autocorrelation(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_rank() {
        let tv = get_test_vec(TestFuncs::Rank);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.rank();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_cumulative_sums_forward() {
        let tv = get_test_vec(TestFuncs::CumulativeSumsForward);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.cumulative_sums_forward();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_cumulative_sums_backward() {
        let tv = get_test_vec(TestFuncs::CumulativeSumsBackward);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.cumulative_sums_backward();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_approximate_entropy() {
        let tv = get_test_vec(TestFuncs::ApproximateEntropy);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.approximate_entropy(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_linear_complexity() {
        let tv = get_test_vec(TestFuncs::LinearComplexity);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.linear_complexity(tv.param);
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_universal() {
        let tv = get_test_vec(TestFuncs::Universal);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.universal();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
    }

    #[test]
    fn test_discrete_fourier() {
        let tv = get_test_vec(TestFuncs::DiscreteFourier);
        let sample: Sample = tv.epsilon_str.into();
        let result = sample.discrete_fourier();
        assert_eq_f64(result.pv, tv.pvalue);
        assert_eq_f64(result.qv, tv.qvalue);
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

        pv.push(sample.serial1(5));
        let now = Instant::now();
        println!("serial: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.serial2(5));
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

        pv.push(sample.longest_run0());
        let now = Instant::now();
        println!("longest_run0: {:.2} s", (now - last).as_secs_f64());
        let last = now;

        pv.push(sample.longest_run1());
        let now = Instant::now();
        println!("longest_run0: {:.2} s", (now - last).as_secs_f64());
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

        pv.push(sample.cumulative_sums_forward());
        let now = Instant::now();
        println!(
            "cumulative_sums_forward: {:.2} s",
            (now - last).as_secs_f64()
        );
        let last = now;

        pv.push(sample.cumulative_sums_backward());
        let now = Instant::now();
        println!(
            "cumulative_sums_backward: {:.2} s",
            (now - last).as_secs_f64()
        );
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
