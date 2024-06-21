use libm::sqrt;
use opengm_rts::*;
use rayon::prelude::*;
use std::fmt::Display;
use std::io::Result;
use std::{env, fs};

const ALPHA: f64 = 0.01;

// GM/T 0005-2021 6.3
const SAMPLE_DISTRIBUTION_K: usize = 10;
const SAMPLE_DISTRIBUTION_ALPHA_T: f64 = 0.0001;

mod tester;
use tester::*;

#[derive(Copy, Clone)]
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

impl Display for TestFuncs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            TestFuncs::Frequency => {
                write!(f, "{:24}", "Frequency")
            }
            TestFuncs::BlockFrequency => {
                write!(f, "{:24}", "BlockFrequency")
            }
            TestFuncs::Poker => {
                write!(f, "{:24}", "Poker")
            }
            TestFuncs::Serial => {
                write!(f, "{:24}", "Serial")
            }
            TestFuncs::Runs => {
                write!(f, "{:24}", "Runs")
            }
            TestFuncs::RunsDistribution => {
                write!(f, "{:24}", "RunsDistribution")
            }
            TestFuncs::LongestRun => {
                write!(f, "{:24}", "LongestRun")
            }
            TestFuncs::BinaryDerivative => {
                write!(f, "{:24}", "BinaryDerivative")
            }
            TestFuncs::Autocorrelation => {
                write!(f, "{:24}", "Autocorrelation")
            }
            TestFuncs::Rank => {
                write!(f, "{:24}", "Rank")
            }
            TestFuncs::CumulativeSums => {
                write!(f, "{:24}", "CumulativeSums")
            }
            TestFuncs::ApproximateEntropy => {
                write!(f, "{:24}", "ApproximateEntropy")
            }
            TestFuncs::LinearComplexity => {
                write!(f, "{:24}", "LinearComplexity")
            }
            TestFuncs::Universal => {
                write!(f, "{:24}", "Universal")
            }
            TestFuncs::DiscreteFourier => {
                write!(f, "{:24}", "DiscreteFourier")
            }
        }
    }
}

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

fn get_testers(funcs: &[TestFuncs], bits: usize) -> Vec<Box<dyn TestFn>> {
    let mut res: Vec<Box<dyn TestFn>> = Vec::new();
    for f in funcs {
        match *f {
            TestFuncs::Frequency => {
                res.push(Box::new(FrequencyTester {}));
            }
            TestFuncs::BlockFrequency => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(BlockFrequencyTester { param }));
                    }
                }
            }
            TestFuncs::Poker => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(PokerTester { param }));
                    }
                }
            }
            TestFuncs::Serial => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(SerialTester { param }));
                    }
                }
            }
            TestFuncs::Runs => {
                res.push(Box::new(RunsTester {}));
            }
            TestFuncs::RunsDistribution => {
                res.push(Box::new(RunsDistributionTester {}));
            }
            TestFuncs::LongestRun => {
                res.push(Box::new(LongestRunTester {}));
            }
            TestFuncs::BinaryDerivative => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(BinaryDerivativeTester { param }));
                    }
                }
            }
            TestFuncs::Autocorrelation => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(AutocorrelationTester { param }));
                    }
                }
            }
            TestFuncs::Rank => {
                res.push(Box::new(RankTester {}));
            }
            TestFuncs::CumulativeSums => res.push(Box::new(CumulativeSumsTester {})),
            TestFuncs::ApproximateEntropy => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(ApproximateEntropyTester { param }));
                    }
                }
            }
            TestFuncs::LinearComplexity => {
                if let Some(params) = get_param(*f, bits) {
                    for param in params {
                        res.push(Box::new(LinearComplexityTester { param }));
                    }
                }
            }
            TestFuncs::Universal => {
                res.push(Box::new(UniversalTester {}));
            }
            TestFuncs::DiscreteFourier => {
                res.push(Box::new(DiscreteFourierTester {}));
            }
        }
    }
    res
}

fn sample_test(s: &Sample, testers: &Vec<Box<dyn TestFn>>) -> Vec<TestResultWrapper> {
    testers.iter().map(|t| t.test_func(s)).collect()
}

fn waterline(alpha: f64, s: usize) -> usize {
    let s = s as f64;
    (s * (1.0 - alpha - 3.0 * sqrt(alpha * (1.0 - alpha) / s))).ceil() as usize
}

fn randomness_test(samples: &Vec<Sample>, test_funcs: &[TestFuncs]) -> (bool, bool) {
    if samples.len() == 0 || test_funcs.len() == 0 {
        println!("No tests to do!");
        return (true, true);
    }
    let testers = get_testers(test_funcs, samples[0].bits());

    let results: Vec<Vec<TestResultWrapper>> = samples
        .par_iter()
        .map(|s| sample_test(s, &testers))
        .collect();
    let line = waterline(ALPHA, samples.len());
    let mut pvalue_res = true;
    let mut qvalue_res: bool = true;

    // results[i][j] is the i-th sample's j'th test
    for j in 0..testers.len() {
        let mut passed = 0;
        let mut q1values = Vec::with_capacity(samples.len());
        let mut q2values = Vec::with_capacity(samples.len());

        for i in 0..samples.len() {
            if results[i][j].result.pass(ALPHA) {
                passed += 1;
            }

            q1values.push(results[i][j].result.qv1);

            if let Some(qv2) = results[i][j].result.qv2 {
                q2values.push(qv2)
            }
        }
        let pt1 = qvalue_distribution(&q1values, SAMPLE_DISTRIBUTION_K);
        let pt2 = qvalue_distribution(&q2values, SAMPLE_DISTRIBUTION_K);

        let this_qvalue_res = pt1 >= SAMPLE_DISTRIBUTION_ALPHA_T && pt2 >= SAMPLE_DISTRIBUTION_ALPHA_T;

        if let Some(param) = results[0][j].param {
            println!(
                "{:2}: Test {} param={:<6}   P-Value: {:>4}/{:<4} {} Q-Value: {}",
                j,
                results[0][j].test_func,
                param,
                passed,
                samples.len(),
                if passed >= line { "Pass,  " } else { "Failed," },
                if q2values.is_empty() {
                    format!("pt={:.6},          {}", pt1, this_qvalue_res)
                } else {
                    format!("pt={:.6} {:.6}, {}", pt1, pt2, this_qvalue_res)
                }
            );
        } else {
            println!(
                "{:2}: Test {}                P-Value: {:>4}/{:<4} {} Q-Value: {}",
                j,
                results[0][j].test_func,
                passed,
                samples.len(),
                if passed >= line { "Pass,  " } else { "Failed," },
                if q2values.is_empty() {
                    format!("pt={:.6},          {}", pt1, this_qvalue_res)
                } else {
                    format!("pt={:.6} {:.6}, {}", pt1, pt2, this_qvalue_res)
                }
            );
        }
        pvalue_res = pvalue_res && (passed >= line);
        qvalue_res = qvalue_res && this_qvalue_res;
    }
    (pvalue_res, qvalue_res)
}

fn read_dir(current_dir: &str) -> Result<Vec<Sample>> {
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

const ALL_TESTS_FUNCS: [TestFuncs; 15] = [
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

const USAGE: &str = r#"opengm_rts v0.1.1
Copyright (c) 2024 The OpenGM Group <opengm@yeah.net>

Usage: ./opengm_rts <path/to/data/dir>

Example: 
$ ls data
000.bin 001.bin 002.bin 003.bin 004.bin ... 999.bin
$ ./opengm_rts ./data
"#;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("{}", USAGE);
        return;
    }
    let result = read_dir(args[1].as_str())
        .expect(&String::from(format!("Read {} error", args[1].as_str())));

    let (p, q) = randomness_test(&result, &ALL_TESTS_FUNCS);
    if p {
        println!("{} P-Value test pass", args[1]);
    } else {
        println!("{} P-Value test failed", args[1]);
    }
    if q {
        println!("{} Q-Value test pass", args[1]);
    } else {
        println!("{} Q-Value test failed", args[1]);
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use rand::RngCore;

    use super::*;

    #[test]
    fn test_rts() {
        let mut samples: Vec<Sample> = Vec::new();
        let mut data = vec![0u8; 1000000 / 8];
        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            rng.fill_bytes(&mut data);
            samples.push(Sample::from(data.as_slice()));
        }
        let start = Instant::now();
        println!("{:?}", randomness_test(&samples, &ALL_TESTS_FUNCS));
        let elapsed = Instant::now() - start;
        println!("Used time: {} s", elapsed.as_secs_f64())
    }

    #[test]
    fn test_waterline() {
        assert_eq!(waterline(ALPHA, 1000), 981);
    }
}
