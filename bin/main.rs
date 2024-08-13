use std::time::Instant;
use opengm_rts::*;
use std::env;

const INFO: &str = r#"opengm_rts v0.1.3
Copyright (c) 2024 The OpenGM Group <opengm@yeah.net>
"#;

const USAGE: &str = r#"Usage: ./opengm_rts <path/to/data/dir>

Example: 
$ ls data
000.bin 001.bin 002.bin 003.bin 004.bin ... 999.bin
$ ./opengm_rts ./data
"#;

fn main() {
    println!("{}", INFO);
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("{}", USAGE);
        return;
    }
    let samples = read_dir(args[1].as_str())
        .expect(&String::from(format!("Read {} error", args[1].as_str())));
    if samples.len() == 0{
        println!("No files in {}!", args[1].as_str());
        return;
    }
    let start_time = Instant::now();

    let n = samples.len();
    let bits = samples[0].bits();
    let w = waterline(ALPHA, samples.len());

    let testers = get_testers(&ALL_TESTS_FUNCS, bits);

    let (presult, qresult) = randomness_test(&samples, &testers);

    let mut pfailed_tests = Vec::new();
    let mut qfailed_tests = Vec::new();

    for tester in testers{
        println!("Test {}:", tester);

        let passed  = count_pvalue_pass(&presult, tester, ALPHA);
        println!("    p_value: {}/{}, {}",passed, n, passed >= w as i32);
        println!("    q_value: {:.4}, {}",qresult.get(&tester).unwrap(),*qresult.get(&tester).unwrap() >= SAMPLE_DISTRIBUTION_ALPHA_T);
        if passed < w as i32{
            pfailed_tests.push(tester);
        }
        if *qresult.get(&tester).unwrap() < SAMPLE_DISTRIBUTION_ALPHA_T{
            qfailed_tests.push(tester);
        }
        println!("");
    }

    if pfailed_tests.len() > 0{
        print!("p_value tests failed: ");
        for t in &pfailed_tests{
            print!("{}", t);
        }
        println!("");
    }
    if qfailed_tests.len() > 0{
        print!("q_value tests failed: ");
        for t in &qfailed_tests{
            print!("{}", t);
        }
        println!("");

    }

    if pfailed_tests.len() == 0 && pfailed_tests.len() == 0{
        println!("Randomness tests passed.")
    }else{
        println!("Randomness tests failed.")
    }
    println!("Randomness tests used time: {:?}",  (Instant::now() - start_time).as_secs_f64());

}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use rand::RngCore;

    use super::*;

    #[test]
    fn test_rts() {
        let testers = get_testers(&ALL_TESTS_FUNCS, 1000000);
        println!("{}", testers);
        let mut samples: Vec<Sample> = Vec::new();
        let mut data = vec![0u8; 1000000 / 8];
        let mut rng = rand::thread_rng();
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
