// Commbined test suites
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Display;
use std::fs;
use std::io::Result;

use crate::Sample;
use crate::*;

pub const ALPHA: f64 = 0.01;
// GM/T 0005-2021 6.3
pub const SAMPLE_DISTRIBUTION_K: usize = 10;
pub const SAMPLE_DISTRIBUTION_ALPHA_T: f64 = 0.0001;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum TestFuncs {
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

impl Display for TestFuncs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            TestFuncs::Frequency => {
                write!(f, "Frequency")
            }
            TestFuncs::BlockFrequency => {
                write!(f, "BlockFrequency")
            }
            TestFuncs::Poker => {
                write!(f, "Poker")
            }
            TestFuncs::Serial => {
                write!(f, "Serial")
            }
            TestFuncs::Runs => {
                write!(f, "Runs")
            }
            TestFuncs::RunsDistribution => {
                write!(f, "RunsDistribution")
            }
            TestFuncs::LongestRun => {
                write!(f, "LongestRun")
            }
            TestFuncs::BinaryDerivative => {
                write!(f, "BinaryDerivative")
            }
            TestFuncs::Autocorrelation => {
                write!(f, "Autocorrelation")
            }
            TestFuncs::Rank => {
                write!(f, "Rank")
            }
            TestFuncs::CumulativeSums => {
                write!(f, "CumulativeSums")
            }
            TestFuncs::ApproximateEntropy => {
                write!(f, "ApproximateEntropy")
            }
            TestFuncs::LinearComplexity => {
                write!(f, "LinearComplexity")
            }
            TestFuncs::Universal => {
                write!(f, "Universal")
            }
            TestFuncs::DiscreteFourier => {
                write!(f, "DiscreteFourier")
            }
        }
    }
}

pub const ALL_TESTS_FUNCS: [TestFuncs; 15] = [
    TestFuncs::Frequency,
    TestFuncs::BlockFrequency,
    TestFuncs::Poker,
    TestFuncs::Serial,
    TestFuncs::Runs,
    TestFuncs::RunsDistribution,
    TestFuncs::LongestRun,
    TestFuncs::BinaryDerivative,
    TestFuncs::Autocorrelation,
    TestFuncs::Rank,
    TestFuncs::CumulativeSums,
    TestFuncs::ApproximateEntropy,
    TestFuncs::LinearComplexity,
    TestFuncs::Universal,
    TestFuncs::DiscreteFourier,
];

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub struct Tester{
    pub f: TestFuncs,
    pub param: Option<i32>,
}

impl Display for Tester{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(param) = self.param{
            write!(f, "{}({})", self.f, param)
        }else{
            write!(f, "{}", self.f)
        }
    }
}

impl Tester {
    fn new(f: TestFuncs, param: Option<i32>) -> Tester {
        Tester{f, param}
    }
    fn test(&self, s: &Sample) -> TestResult{
        let mut m = 0;
        if self.param.is_some(){
            m = self.param.unwrap();
        }

        match self.f{
            TestFuncs::Frequency => s.frequency(),
            TestFuncs::BlockFrequency => s.block_frequency(m),
            TestFuncs::Poker => s.poker(m),
            TestFuncs::Serial => s.serial(m),
            TestFuncs::Runs => s.runs(),
            TestFuncs::RunsDistribution => s.runs_distribution(),
            TestFuncs::LongestRun => s.longest_run(),
            TestFuncs::BinaryDerivative => s.binary_derivative(m),
            TestFuncs::Autocorrelation => s.autocorrelation(m),
            TestFuncs::Rank => s.rank(),
            TestFuncs::CumulativeSums => s.cumulative_sums(),
            TestFuncs::ApproximateEntropy => s.approximate_entropy(m),
            TestFuncs::LinearComplexity => s.linear_complexity(m),
            TestFuncs::Universal => s.universal(),
            TestFuncs::DiscreteFourier => s.discrete_fourier(),
        }
    }
}

// #[derive(Debug)]
// pub struct TestResultWrapper {
//     pub result: TestResult,
//     pub test_func: TestFuncs,
//     pub param: Option<i32>,
// }


fn get_param(f: TestFuncs, bits: usize) -> Option<Vec<i32>> {
    if bits <= 20000 {
        match f {
            TestFuncs::BlockFrequency => Some(vec![1000]),
            TestFuncs::Poker => Some(vec![4, 8]),
            TestFuncs::Serial => Some(vec![3, 5]),
            TestFuncs::LongestRun => Some(vec![128]),
            TestFuncs::BinaryDerivative => Some(vec![3, 7]),
            TestFuncs::Autocorrelation => Some(vec![2, 8, 16]),
            TestFuncs::ApproximateEntropy => Some(vec![2, 5]),
            _ => None,
        }
    } else if bits <= 1000000 {
        // Universal has fixed parameters L,Q
        match f {
            TestFuncs::BlockFrequency => Some(vec![10000]),
            TestFuncs::Poker => Some(vec![4, 8]),
            TestFuncs::Serial => Some(vec![3, 5]),
            TestFuncs::LongestRun => Some(vec![10000]),
            TestFuncs::BinaryDerivative => Some(vec![3, 7]),
            TestFuncs::Autocorrelation => Some(vec![1, 2, 8, 16]),
            TestFuncs::ApproximateEntropy => Some(vec![2, 5]),
            TestFuncs::LinearComplexity => Some(vec![500, 1000]),
            _ => None,
        }
    } else {
        match f {
            TestFuncs::BlockFrequency => Some(vec![1000000]),
            TestFuncs::Poker => Some(vec![4, 8]),
            TestFuncs::Serial => Some(vec![3, 5, 7]),
            TestFuncs::LongestRun => Some(vec![10000]),
            TestFuncs::BinaryDerivative => Some(vec![3, 7, 15]),
            TestFuncs::Autocorrelation => Some(vec![1, 2, 8, 16, 32]),
            TestFuncs::ApproximateEntropy => Some(vec![2, 5]),
            TestFuncs::LinearComplexity => Some(vec![5000]),
            _ => None,
        }
    }
}

pub fn get_testers(funcs: &[TestFuncs], bits: usize) -> Vec<Tester> {
    let mut res: Vec<_> = Vec::new();
    for f in funcs {
        if let Some(params) = get_param(*f, bits) {
            for param in params {
                res.push(Tester {f:*f, param:Some(param)});
            }
        }else{
            res.push(Tester {f:*f, param:None});
        }
    }
    res
}

pub fn waterline(alpha: f64, s: usize) -> usize {
    let s = s as f64;
    (s * (1.0 - alpha - 3.0 * (alpha * (1.0 - alpha) / s).sqrt())).ceil() as usize
}

fn sample_test(
    s: &Sample,
    testers: &[Tester],
) -> HashMap<Tester, TestResult> {
    let mut result = HashMap::<Tester, TestResult>::new();
    for t in testers {
        result.insert(*t, t.test(s));
    }
    return result;
}

// return p_value result and q_value result.
// p_value result returns [samples[0].Result, samples[1].Result,...]
// samples[i].Result: Tester: TestResult.
// q_value result returns hashmap: Tester: q_value distribution.
pub fn randomness_test(
    samples: &Vec<Sample>,
    testers: &[Tester],
) -> (
    Vec<HashMap<Tester, TestResult>>,
    HashMap<Tester, f64>,
) {
    if samples.len() == 0 || testers.len() == 0 {
        return (
            Vec::<HashMap<Tester, TestResult>>::new(),
            HashMap::new(),
        );
    }
    let presult = samples
        .par_iter() // multiple threads by rayon
        // .iter() // sigle thread
        .map(|s| sample_test(s, &testers))
        .collect::<Vec<HashMap<Tester, TestResult>>>();

    let mut qresult = HashMap::with_capacity(testers.len());
    for f in testers{
        qresult.insert(*f, compute_qvalue_distribution(&presult, *f));
    }
    (presult, qresult)
}

// let line = waterline(ALPHA, samples.len());
pub fn count_pvalue_pass(
    res: &Vec<HashMap<Tester, TestResult>>,
    f: Tester,
    alpha: f64,
) -> i32 {
    res.iter()
        .map(|t| if t[&f].pass(alpha) { 1 } else { 0 })
        .sum()
}

fn compute_qvalue_distribution(
    res: &Vec<HashMap<Tester, TestResult>>,
    f: Tester,
) -> f64 {
    let mut qvalues = Vec::with_capacity(res.len() * 2);

    for r in res {
        let result = &r.get(&f).unwrap();
        qvalues.push(result.qv1);
        if let Some(qv2) = result.qv2 {
            qvalues.push(qv2)
        }
    }
    qvalue_distribution(&qvalues, SAMPLE_DISTRIBUTION_K)
}

pub fn read_dir(current_dir: &str) -> Result<Vec<Sample>> {
    let mut result = Vec::new();
    for entry in fs::read_dir(current_dir)? {
        let entry = entry?;
        let path = entry.path();

        let metadata = fs::metadata(&path)?;
        if metadata.is_file() {
            if let Ok(content) = fs::read(path) {
                result.push(Sample::from(content.as_slice()));
            }
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use rand::RngCore;

    use super::*;

    // cargo test --release --package opengm_rts --lib -- tester::tests::test_rts --exact --show-output
    #[test]
    fn test_rts() {
        let mut samples: Vec<Sample> = Vec::new();
        let bits = 1000000;
        let mut data = vec![0u8; bits/8];
        let mut rng = rand::thread_rng();
        let testers = get_testers(&ALL_TESTS_FUNCS, bits);

        for _ in 0..1000 {
            rng.fill_bytes(&mut data);
            samples.push(Sample::from(data.as_slice()));
        }
        let start = Instant::now();
        println!("{:?}", randomness_test(&samples, &testers));
        let elapsed = Instant::now() - start;
        println!("Used time: {} s", elapsed.as_secs_f64())
    }

    #[test]
    fn test_waterline() {
        assert_eq!(waterline(ALPHA, 1000), 981);
    }
}