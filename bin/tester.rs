use opengm_rts::*;

use crate::TestFuncs;

pub struct TestResultWrapper {
   pub result: TestResult,
   pub test_func: TestFuncs,
   pub param: Option<i32>,
}

pub trait TestFn: Sync+Send {
    fn test_func(&self, s: &Sample) -> TestResultWrapper;
}

pub struct FrequencyTester {}

impl TestFn for FrequencyTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.frequency(),
            test_func: TestFuncs::Frequency,
            param: None,
        }
    }
}

pub struct BlockFrequencyTester {
    pub param: i32,
}

impl TestFn for BlockFrequencyTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.block_frequency(self.param),
            test_func: TestFuncs::BlockFrequency,
            param: Some(self.param),
        }
    }
}

pub struct PokerTester {
    pub param: i32,
}
impl TestFn for PokerTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.poker(self.param),
            test_func: TestFuncs::Poker,
            param: Some(self.param),
        }
        
    }
}

pub struct SerialTester {
    pub param: i32,
}
impl TestFn for SerialTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.serial(self.param),
            test_func: TestFuncs::Serial,
            param: Some(self.param),
        }
    }
}

pub struct RunsTester {}
impl TestFn for RunsTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.runs(),
            test_func: TestFuncs::Runs,
            param: None,
        }
    }
}

pub struct RunsDistributionTester {}
impl TestFn for RunsDistributionTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.runs_distribution(),
            test_func: TestFuncs::RunsDistribution,
            param: None,
        }
    }
}

pub struct LongestRunTester {}
impl TestFn for LongestRunTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.longest_run(),
            test_func: TestFuncs::LongestRun,
            param: None,
        }
    }
}

pub struct BinaryDerivativeTester {
    pub param: i32,
}
impl TestFn for BinaryDerivativeTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.binary_derivative(self.param),
            test_func: TestFuncs::BinaryDerivative,
            param: Some(self.param),
        }
        
    }
}

pub struct AutocorrelationTester {
    pub param: i32,
}
impl TestFn for AutocorrelationTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.autocorrelation(self.param),
            test_func: TestFuncs::Autocorrelation,
            param: Some(self.param),
        }
        
    }
}

pub struct RankTester {}
impl TestFn for RankTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.rank(),
            test_func: TestFuncs::Rank,
            param: None,
        }
        
    }
}

pub struct CumulativeSumsTester {}
impl TestFn for CumulativeSumsTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.cumulative_sums(),
            test_func: TestFuncs::CumulativeSums,
            param: None,
        }
        
    }
}

pub struct ApproximateEntropyTester {
    pub param: i32,
}
impl TestFn for ApproximateEntropyTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.approximate_entropy(self.param),
            test_func: TestFuncs::ApproximateEntropy,
            param: Some(self.param),
        }
        
    }
}

pub struct LinearComplexityTester {
    pub param: i32,
}
impl TestFn for LinearComplexityTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.linear_complexity(self.param),
            test_func: TestFuncs::LinearComplexity,
            param: Some(self.param),
        }
        
    }
}

pub struct UniversalTester {}
impl TestFn for UniversalTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.universal(),
            test_func: TestFuncs::Universal,
            param: None,
        }
        
    }
}

pub struct DiscreteFourierTester {}
impl TestFn for DiscreteFourierTester {
    fn test_func(&self, s: &Sample) -> TestResultWrapper {
        TestResultWrapper {
            result: s.discrete_fourier(),
            test_func: TestFuncs::DiscreteFourier,
            param: None,
        }
    }
}

