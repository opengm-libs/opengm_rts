#![allow(dead_code)]
#![allow(non_snake_case)]
#![feature(test)]

mod tester;
mod testsuits;
use std::sync::Mutex;

pub use tester::*;
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

impl PartialEq for TestResult {
    fn eq(&self, other: &Self) -> bool {
        (self.pv - other.pv).abs() < 0.0001
            && (self.qv - other.qv).abs() < 0.0001
    }
}

impl TestResult {
    pub fn pass(&self, alpha: f64) -> bool {
        self.pv >= alpha
    }
}

const MAX_OVERLAPPING_PATTERN: usize = 7;

pub struct Sample {
    // Note: for a byte x, bit 0 is (x>>7) & 1 and bit 7 is x & 1.
    // Note: for a u64 x, bit 0 is (x>>63) & 1 and bit 63 is x & 1.
    b: Vec<u8>,
    b64: Vec<u64>,
    bit_length: usize,

    #[cfg(test)]
    e: Vec<u8>,

    // popcount for e
    pop: u64,

    // 重叠子序列和近似熵
    patterns: Mutex<Vec<Option<Vec<u64>>>>,
    max_patterns: usize,
}

fn combine(b: &[u8]) -> Vec<u64> {
    let mut b64 = Vec::with_capacity((b.len() + 7) / 8);
    let full_chunks = b.len() & (!7);
    for chunk in b[..full_chunks].chunks_exact(8) {
        b64.push(u64::from_be_bytes(chunk.try_into().unwrap()));
    }

    if full_chunks < b.len() {
        b64.push(u64_from_be_slice(&b[full_chunks..]));
    }
    b64
}

// from bit string("10100...")
impl From<&str> for Sample {
    fn from(bit_string: &str) -> Self {
        let mut e = Vec::with_capacity(bit_string.len());
        let mut b = Vec::with_capacity((bit_string.len() + 7) / 8);
        for c in bit_string.chars() {
            e.push(c as u8 - '0' as u8);
        }

        let full_chunks = e.len() & (!7);
        for chunk in e[..full_chunks].chunks_exact(8) {
            let mut x = 0;
            for i in chunk {
                x = (x << 1) | ((*i) as u8);
            }
            b.push(x);
        }

        if full_chunks < e.len() {
            let mut c = 0;
            for (i, x) in e[full_chunks..].iter().enumerate() {
                c |= (*x) << (7 - i);
            }
            b.push(c);
        }
        let b64 = combine(&b);
        let pop = popcount(&b);
        Sample {
            b: b,
            b64: b64,
            #[cfg(test)]
            e: e,
            bit_length: bit_string.len(),
            pop: pop,
            patterns: Mutex::new(vec![None; MAX_OVERLAPPING_PATTERN+1]),
            max_patterns: MAX_OVERLAPPING_PATTERN,
        }
    }
}

impl From<&[u8]> for Sample {
    fn from(byte_slice: &[u8]) -> Self {
        #[cfg(test)]
        let mut e = Vec::with_capacity(byte_slice.len() * 8);
        #[cfg(test)]
        for c in byte_slice {
            for j in 0..8 {
                e.push((c >> (7 - j)) & 1)
            }
        }

        let b64 = combine(&byte_slice);
        let pop = popcount(&byte_slice);
        Sample {
            b: byte_slice.to_vec(),
            b64,
            bit_length: byte_slice.len() * 8,
            #[cfg(test)]
            e,
            pop: pop,
            patterns:Mutex::new(vec![None; MAX_OVERLAPPING_PATTERN+1]),
            max_patterns: MAX_OVERLAPPING_PATTERN,
        }
    }
}

#[inline(always)]
fn get_bit_unchecked(b64: &[u64], i: usize) -> u8 {
    unsafe { (b64.get_unchecked(i / 64) >> (63 - (i & 63))) as u8 & 1 }
}

impl Sample {
    pub fn new(bytes: &[u8]) -> Sample {
        return bytes.into();
    }

    #[inline(always)]
    pub fn get_bit_unchecked(&self, i: usize) -> u8 {
        get_bit_unchecked(&self.b64, i)
    }

    /// returns the number of bits of sample.
    pub fn len(&self) -> usize {
        self.bit_length
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
pub mod tests {
    use std::collections::HashMap;
    use std::time::Instant;
    use rand::RngCore;
    use super::test_data::E;
    use super::*;


    // Test for all 15 test functions.
    // cargo test --release --package opengm_rts --lib -- tests::test_all --exact --show-output
    #[test]
    fn test_all() {
        let sample: Sample = E.into();
        let mut pv = HashMap::new();

        let start = Instant::now();
        let last = start;

        pv.insert(TestFuncs::Frequency, sample.frequency());
        let now = Instant::now();
        println!("frequency: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::BlockFrequency, sample.block_frequency(10000));
        let now = Instant::now();
        println!("block_frequency: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Poker, sample.poker(8));
        let now = Instant::now();
        println!("poker: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Serial1, sample.serial1(5));
        let now = Instant::now();
        println!("serial: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Serial2, sample.serial2(5));
        let now = Instant::now();
        println!("serial: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Runs, sample.runs());
        let now = Instant::now();
        println!("runs: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::RunsDistribution, sample.runs_distribution());
        let now = Instant::now();
        println!("runs_distribution: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::LongestRun0, sample.longest_run0());
        let now = Instant::now();
        println!("longest_run0: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::LongestRun1, sample.longest_run1());
        let now = Instant::now();
        println!("longest_run0: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::BinaryDerivative, sample.binary_derivative(7));
        let now = Instant::now();
        println!("binary_derivative: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Autocorrelation, sample.autocorrelation(16));
        let now = Instant::now();
        println!("autocorrelation: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Rank, sample.rank());
        let now = Instant::now();
        println!("rank: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(
            TestFuncs::CumulativeSumsForward,
            sample.cumulative_sums_forward(),
        );
        let now = Instant::now();
        println!(
            "cumulative_sums_forward: {:.6} s",
            (now - last).as_secs_f64()
        );
        let last = now;

        pv.insert(
            TestFuncs::CumulativeSumsBackward,
            sample.cumulative_sums_backward(),
        );
        let now = Instant::now();
        println!(
            "cumulative_sums_backward: {:.6} s",
            (now - last).as_secs_f64()
        );
        let last = now;

        pv.insert(TestFuncs::ApproximateEntropy, sample.approximate_entropy(5));
        let now = Instant::now();
        println!("approximate_entropy: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::LinearComplexity, sample.linear_complexity(1000));
        let now = Instant::now();
        println!("linear_complexity: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::Universal, sample.universal());
        let now = Instant::now();
        println!("universal: {:.6} s", (now - last).as_secs_f64());
        let last = now;

        pv.insert(TestFuncs::DiscreteFourier, sample.discrete_fourier());
        let now = Instant::now();
        println!("discrete_fourier: {:.6} s", (now - last).as_secs_f64());

        println!("total: {:.6} s", (Instant::now() - start).as_secs_f64());
    }

    // Test for all 15 test functions.
    // cargo test --release --package opengm_rts --lib -- tests::test_all_time --exact --show-output
    #[test]
    fn test_all_time() {
        let mut samples: Vec<Sample> = Vec::new();
        let mut data = vec![0u8; 1000000/8];
        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            rng.fill_bytes(&mut data);
            samples.push(Sample::from(data.as_slice()));
        }
        let mut pv = HashMap::new();

        // total: 74.3s
        let start = Instant::now();
        for sample in samples {
            // 0s
            // pv.insert(TestFuncs::Frequency, sample.frequency());

            // // 0.01
            // pv.insert(TestFuncs::BlockFrequency, sample.block_frequency(10000));

            // pv.insert(TestFuncs::Poker, sample.poker(8));

            // 0.90 s
            // pv.insert(TestFuncs::Serial1, sample.serial1(7));

            // 1.47 s
            // pv.insert(TestFuncs::Serial2, sample.serial2(7));

            // 0.02 s
            // pv.insert(TestFuncs::Runs, sample.runs());

            // 4.21 s - TODO
            // pv.insert(TestFuncs::RunsDistribution, sample.runs_distribution());

            // 0.96s - TODO
            // pv.insert(TestFuncs::LongestRun0, sample.longest_run0());
            // pv.insert(TestFuncs::LongestRun1, sample.longest_run1());

            // 0.83s, 0.04 s
            // pv.insert(TestFuncs::BinaryDerivative, sample.binary_derivative(7));

            // 0.65s, 0.02s
            // pv.insert(TestFuncs::Autocorrelation, sample.autocorrelation(8));

            // 6.67s -> 1.89 s
            // pv.insert(TestFuncs::Rank, sample.rank());

            // 1s
            // pv.insert(
            //     TestFuncs::CumulativeSumsForward,
            //     sample.cumulative_sums_forward(),
            // );

            // 1s
            // pv.insert(
            //     TestFuncs::CumulativeSumsBackward,
            //     sample.cumulative_sums_backward(),
            // );

            // 4.69s, 0.94s
            // pv.insert(
            //     TestFuncs::ApproximateEntropy,
            //     sample.approximate_entropy(5),
            // );

            // m=1000:20s
            // m=500:15s
            // pv.insert(
            //     TestFuncs::LinearComplexity,
            //     sample.linear_complexity(500),
            // );

            // 0.58s
            // pv.insert(TestFuncs::Universal, sample.universal());

            // 27s
            pv.insert(TestFuncs::DiscreteFourier, sample.discrete_fourier());
        }
        println!("total: {:.2} s", (Instant::now() - start).as_secs_f64());
    }
}
