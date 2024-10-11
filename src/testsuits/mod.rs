pub(crate) mod util;

mod approximate_entropy;
mod autocorrelation;
mod binary_derivative;
mod block_frequency;
mod cumulative_sums;
mod discrete_fourier;
mod frequency;
mod linear_complexity;
pub(crate) mod linear_complexity_simd;
mod longest_run;
mod poker;
mod rank;
mod rank_u64;
mod runs;
mod runs_distribution;
mod serial;
mod universal;

pub(crate) use approximate_entropy::*;
pub(crate) use autocorrelation::*;
pub(crate) use binary_derivative::*;
pub(crate) use block_frequency::*;
pub(crate) use cumulative_sums::*;
pub(crate) use discrete_fourier::*;
pub(crate) use frequency::*;
pub(crate) use linear_complexity::*;
pub(crate) use longest_run::*;
pub(crate) use poker::*;
pub(crate) use rank::*;
pub(crate) use runs::*;
pub(crate) use runs_distribution::*;
pub(crate) use serial::*;
pub(crate) use universal::*;



//Test vectors
#[cfg(test)]
mod tests {
    use crate::{test_data::E, Sample, TestFuncs, TestResult};

    pub fn get_test_vec(f: TestFuncs) -> (&'static str, i32, TestResult) {
        match f{
        // 1.frequency 频度检测
        TestFuncs::Frequency => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            0,
            TestResult{
            pv: 0.215925,
            qv: 0.892038,
            }
        ),
        // 2.block_frequency 块内频度检测
        TestFuncs::BlockFrequency => (
            "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            10,
            TestResult{
            pv : 0.706438,
            qv: 0.706438,
            }
        ),

        // 3.poker 扑克检测
        TestFuncs::Poker => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            4,
            TestResult{
            pv : 0.213734,
            qv: 0.213734,
            }
        ),
        // 4.serial
        TestFuncs::Serial1 => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            2,
            TestResult{
            pv : 0.436868,
            qv:  0.436868,
            }
        ),
        // 4.serial
        TestFuncs::Serial2 => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            2,
            TestResult{
            pv : 0.723674,
            qv: 0.723674,
            }
        ),
        // 5.runs 游程总数检测
        TestFuncs::Runs => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            0,
            TestResult{
            pv : 0.620729,
            qv: 0.310364,
            }
        ),
        // 6)runs_distribution 游程分布检测
        TestFuncs::RunsDistribution => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            0,
            TestResult{
            pv : 0.970152,
            qv: 0.970152,
            }
        ),
        // 7.longest_run 块内最大0游程检测
        TestFuncs::LongestRun0 => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            0,
            TestResult{
            pv : 0.839299,
            qv: 0.839299,
            },
        ),
        // 7)longest_run 块内最大1游程检测
        TestFuncs::LongestRun1 => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            0,
            TestResult{
            pv : 0.180598,
            qv: 0.180598,
            }
        ),
        // 8.binary_derivative 二元推导检测
        TestFuncs::BinaryDerivative => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            3,
            TestResult{
            pv : 0.039669,
            qv: 0.980166,
        },
    ),
        // 9)auto_correlation 自相关检测
        TestFuncs::Autocorrelation => (
            "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010",
            1,
            TestResult{
            pv : 0.790080,
            qv: 0.395040,
            }
        ),
        // 10.rank 矩阵秩检测
        TestFuncs::Rank => (
            E,
            0,
            TestResult{
            pv : 0.307543,
            qv: 0.307543,
            }
        ),
        // 11.cumulative_sums 前向累加和检测
        TestFuncs::CumulativeSumsForward => (
            "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            0,
            TestResult{
            pv : 0.219194,
            qv:  0.219194,
            }
        ),
        // 11.cumulative_sums 后向累加和检测
        TestFuncs::CumulativeSumsBackward => (
            "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            0,
            TestResult{
            pv : 0.114866,
            qv:  0.114866,
            }
        ),
        // 12.approximate_entropy 近似熵检测
        TestFuncs::ApproximateEntropy => (
            "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            2,
            TestResult{
            pv : 0.235301,
            qv: 0.235301,
            }
        ),
        // 13.linear_complexity 线性复杂度检测
        TestFuncs::LinearComplexity => (
            E,
            1000,
            TestResult{
            pv : 0.844721,
            qv: 0.844721,
            }
        ),
        // 14.universal Maurer通用统计检测
        TestFuncs::Universal => (
            E,
            0,
            TestResult{
            pv : 0.282568,
            qv: 0.141284,
            }
        ),
        // 15.discrete_fourier_transform 离散傅立叶检测
        TestFuncs::DiscreteFourier => (
            "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000",
            0,
            TestResult{
            pv: 0.654721,
            qv: 0.327360,
            }
        ),
    }
    }

    pub fn get_test_vec_e(f: TestFuncs) -> (i32, TestResult) {
        match f {
            TestFuncs::Frequency => (
                0,
                TestResult {
                    pv: 0.9537486285283232,
                    qv: 0.4768743142641616,
                },
            ),

            TestFuncs::BlockFrequency => (
                10000,
                TestResult {
                    pv: 0.6762267238527553,
                    qv: 0.6762267238527553,
                },
            ),
            TestFuncs::Poker => (
                8,
                TestResult {
                    pv: 0.023946720623714876,
                    qv: 0.023946720623714876,
                },
            ),

            TestFuncs::Serial1 => (
                5,
                TestResult {
                    pv: 0.22578270155911137,
                    qv: 0.22578270155911137,
                },
            ),
            TestFuncs::Serial2 => (
                5,
                TestResult {
                    pv: 0.05749920567875454,
                    qv: 0.05749920567875454,
                },
            ),

            TestFuncs::Runs => (
                0,
                TestResult {
                    pv: 0.5619168850302544,
                    qv: 0.7190415574848728,
                },
            ),

            TestFuncs::RunsDistribution => (
                0,
                TestResult {
                    pv: 0.7724118883505677,
                    qv: 0.7724118883505677,
                },
            ),

            TestFuncs::LongestRun0 => (
                0,
                TestResult {
                    pv: 0.43786054706891864,
                    qv: 0.43786054706891864,
                },
            ),

            TestFuncs::LongestRun1 => (
                0,
                TestResult {
                    pv: 0.7183547821378014,
                    qv: 0.7183547821378014,
                },
            ),

            TestFuncs::BinaryDerivative => (
                7,
                TestResult {
                    pv: 0.7603653392831897,
                    qv: 0.6198173303584051,
                },
            ),

            TestFuncs::Autocorrelation => (
                16,
                TestResult {
                    pv: 0.912408677003725,
                    qv: 0.5437956614981375,
                },
            ),

            TestFuncs::Rank => (
                0,
                TestResult {
                    pv: 0.30754349438978223,
                    qv: 0.30754349438978223,
                },
            ),

            TestFuncs::CumulativeSumsBackward => (
                0,
                TestResult {
                    pv: 0.7242653099698069,
                    qv: 0.7242653099698069,
                },
            ),

            TestFuncs::CumulativeSumsForward => (
                0,
                TestResult {
                    pv: 0.6698864641681426,
                    qv: 0.6698864641681426,
                },
            ),

            TestFuncs::ApproximateEntropy => (
                5,
                TestResult {
                    pv: 0.3616876318995624,
                    qv: 0.3616876318995624,
                },
            ),

            TestFuncs::LinearComplexity => (
                1000,
                TestResult {
                    pv: 0.8447206463007337,
                    qv: 0.8447206463007337,
                },
            ),

            TestFuncs::Universal => (
                0,
                TestResult {
                    pv: 0.28256794782574396,
                    qv: 0.14128397391287198,
                },
            ),
            TestFuncs::DiscreteFourier => (
                0,
                TestResult {
                    pv: 0.8510101449582823,
                    qv: 0.42550507247914116,
                },
            ),
        }
    }

    #[test]
    fn test_frequency() {
        let tv = get_test_vec(TestFuncs::Frequency);
        let sample: Sample = tv.0.into();
        let result = sample.frequency();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_block_frequency() {
        let tv = get_test_vec(TestFuncs::BlockFrequency);
        let sample: Sample = tv.0.into();
        assert_eq!(sample.block_frequency(tv.1), tv.2);
    }

    #[test]
    fn test_poker() {
        let tv = get_test_vec(TestFuncs::Poker);
        let sample: Sample = tv.0.into();
        let result = sample.poker(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_serial1() {
        let tv = get_test_vec(TestFuncs::Serial1);
        let sample: Sample = tv.0.into();
        let result = sample.serial1(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_serial2() {
        let tv = get_test_vec(TestFuncs::Serial2);
        let sample: Sample = tv.0.into();
        let result = sample.serial2(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_runs() {
        let tv = get_test_vec(TestFuncs::Runs);
        let sample: Sample = tv.0.into();
        let result = sample.runs();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_runs_distribution() {
        let tv = get_test_vec(TestFuncs::RunsDistribution);
        let sample: Sample = tv.0.into();
        let result = sample.runs_distribution();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_longest_run0() {
        let tv = get_test_vec(TestFuncs::LongestRun0);
        let sample: Sample = tv.0.into();
        let result = sample.longest_run0();
        assert_eq!(result, tv.2);
    }
    #[test]
    fn test_longest_run1() {
        let tv = get_test_vec(TestFuncs::LongestRun1);
        let sample: Sample = tv.0.into();
        let result = sample.longest_run1();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_binary_derivative() {
        let tv = get_test_vec(TestFuncs::BinaryDerivative);
        let sample: Sample = tv.0.into();
        let result = sample.binary_derivative(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_autocorrelation() {
        let tv = get_test_vec(TestFuncs::Autocorrelation);
        let sample: Sample = tv.0.into();
        let result = sample.autocorrelation(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_rank() {
        let tv = get_test_vec(TestFuncs::Rank);
        let sample: Sample = tv.0.into();
        let result = sample.rank();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_cumulative_sums_forward() {
        let tv = get_test_vec(TestFuncs::CumulativeSumsForward);
        let sample: Sample = tv.0.into();
        let result = sample.cumulative_sums_forward();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_cumulative_sums_backward() {
        let tv = get_test_vec(TestFuncs::CumulativeSumsBackward);
        let sample: Sample = tv.0.into();
        let result = sample.cumulative_sums_backward();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_approximate_entropy() {
        let tv = get_test_vec(TestFuncs::ApproximateEntropy);
        let sample: Sample = tv.0.into();
        let result = sample.approximate_entropy(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_linear_complexity() {
        let tv = get_test_vec(TestFuncs::LinearComplexity);
        let sample: Sample = tv.0.into();
        let result = sample.linear_complexity(tv.1);
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_universal() {
        let tv = get_test_vec(TestFuncs::Universal);
        let sample: Sample = tv.0.into();
        let result = sample.universal();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_discrete_fourier() {
        let tv = get_test_vec(TestFuncs::DiscreteFourier);
        let sample: Sample = tv.0.into();
        let result = sample.discrete_fourier();
        assert_eq!(result, tv.2);
    }

    #[test]
    fn test_all_fn() {
        let sample: Sample = E.into();
        assert_eq!(get_test_vec_e(TestFuncs::Frequency).1, sample.frequency());
        assert_eq!(
            get_test_vec_e(TestFuncs::BlockFrequency).1,
            sample.block_frequency(10000)
        );
        assert_eq!(get_test_vec_e(TestFuncs::Poker).1, sample.poker(8));
        assert_eq!(get_test_vec_e(TestFuncs::Serial1).1, sample.serial1(5));
        assert_eq!(get_test_vec_e(TestFuncs::Serial2).1, sample.serial2(5));
        assert_eq!(get_test_vec_e(TestFuncs::Runs).1, sample.runs());
        assert_eq!(
            get_test_vec_e(TestFuncs::RunsDistribution).1,
            sample.runs_distribution()
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::LongestRun0).1,
            sample.longest_run0()
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::LongestRun1).1,
            sample.longest_run1()
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::BinaryDerivative).1,
            sample.binary_derivative(7)
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::Autocorrelation).1,
            sample.autocorrelation(16)
        );
        assert_eq!(get_test_vec_e(TestFuncs::Rank).1, sample.rank());
        assert_eq!(
            get_test_vec_e(TestFuncs::CumulativeSumsForward).1,
            sample.cumulative_sums_forward()
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::CumulativeSumsBackward).1,
            sample.cumulative_sums_backward()
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::ApproximateEntropy).1,
            sample.approximate_entropy(5)
        );
        assert_eq!(
            get_test_vec_e(TestFuncs::LinearComplexity).1,
            sample.linear_complexity(1000)
        );
        assert_eq!(get_test_vec_e(TestFuncs::Universal).1, sample.universal());
        assert_eq!(
            get_test_vec_e(TestFuncs::DiscreteFourier).1,
            sample.discrete_fourier()
        );
    }
}
