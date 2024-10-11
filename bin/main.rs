use std::time::Instant;
use opengm_rts::*;
use std::env;

const INFO: &str = r#"opengm_rts v0.2.2
Copyright (c) 2024 The OpenGM Group <opengm@yeah.net>
"#;

const USAGE: &str = r#"Usage: ./opengm_rts <path/to/data/dir>

Example: 
$ ls data
000.bin 001.bin 002.bin 003.bin 004.bin ... 999.bin
$ ./opengm_rts ./data
"#;

macro_rules! pprint {
    ($v: expr, $n: expr) => {
        {
            let s = format!("{}", $v);
            assert!(s.len() <= $n);
            print!("{}", s);
            (0..($n-s.len())).for_each(|_| print!(" "));
        }
    };
}


macro_rules! pprint_center {
    ($v: expr, $n: expr) => {
        {
            let s = format!("{}", $v);
            assert!(s.len() <= $n);
            (0..($n-s.len())/2).for_each(|_| print!(" "));
            print!("{}", s);
            (0..($n-s.len() - ($n-s.len())/2)).for_each(|_| print!(" "));
        }
    };
}

const COL1 :usize = 24;
const COL2: usize = 12;
const COL3: usize = 6;
const COL4: usize = 12;
const COL5: usize = 6;
const COL: usize = COL1 + COL2 + COL3 + COL4 + COL5 + 4;

const PASS: &str = "PASS";
const FAIL: &str = "FAIL";


fn print_line(delimiter: &str){
    print!("+");
    (0..COL1).for_each(|_| print!("-"));
    print!("{}", delimiter);
    (0..COL2).for_each(|_| print!("-"));
    print!("{}", delimiter);
    (0..COL3).for_each(|_| print!("-"));
    print!("{}", delimiter);
    (0..COL4).for_each(|_| print!("-"));
    print!("{}", delimiter);
    (0..COL5).for_each(|_| print!("-"));
    println!("+");
}
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
    let bits = samples[0].len();
    let w = waterline(ALPHA, samples.len());

    let testers = get_testers(&ALL_TESTS_FUNCS, bits);

    let (presult, qresult) = randomness_test(&samples, &testers);

    let mut pfailed_tests = Vec::new();
    let mut qfailed_tests = Vec::new();

    print_line("-");
    print!("|");pprint!(format!("Number of samples: {}", n), COL);println!("|");
    print!("|");pprint!(format!("Bits per sample:   {}", bits), COL);println!("|");
    print!("|");pprint!(format!("P_value threshold: {}", w), COL);println!("|");
    print!("|");pprint!(format!("Q_value threshold: {}", SAMPLE_DISTRIBUTION_ALPHA_T), COL);println!("|");
    
    print_line("-");

    print!("|");
    pprint!(" ", COL1); print!(" ");
    pprint_center!("p_value", COL2); print!(" ");
    pprint!(" ", COL3); print!(" ");
    pprint_center!("q_value", COL4); print!(" ");
    pprint!(" ", COL5);println!("|");

    print_line("+");

    for tester in testers{
        print!("|");
        pprint!(tester, COL1);
        print!("|");

        let passed  = count_pvalue_pass(&presult, tester, ALPHA);
        pprint_center!(format!("{}/{}", passed, n), COL2);
        print!("|");

        let pv_result = passed >= w as i32;
        pprint_center!(if pv_result {PASS} else {FAIL}, COL3);
        print!("|");

        pprint_center!(format!("{:.4}",qresult.get(&tester).unwrap()), COL4);
        print!("|");
        let qv_result = *qresult.get(&tester).unwrap() >= SAMPLE_DISTRIBUTION_ALPHA_T;
        pprint_center!(if qv_result {PASS} else {FAIL}, COL5);
        println!("|");

        if !pv_result {
            pfailed_tests.push(tester);
        }
        if !qv_result{
            qfailed_tests.push(tester);
        }
    }

    print_line("+");
    
    // if pfailed_tests.len() > 0{
    //     print!("p_value测试失败的检测有: {:?}\n\n", pfailed_tests.iter().map(|t| format!("{}",t)).collect::<Vec<_>>());
    //     print_line("+");
    // }


    // if qfailed_tests.len() > 0{
    //     print!("q_value测试失败的检测有: {:?}\n\n", qfailed_tests.iter().map(|t| format!("{}",t)).collect::<Vec<_>>());
    //     print_line("+");
    // }


    if pfailed_tests.len() == 0 && pfailed_tests.len() == 0{
        print!("|");
        pprint!("Randomness test PASS.", COL);
        println!("|");
    }else{
        print!("|");
        println!("Randomness test FAIL.");
        println!("|");
    }

    print!("|");
    pprint!(format!("Used time: {:.01} seconds",  (Instant::now() - start_time).as_secs_f64()), COL);
    println!("|");
    print_line("-");
}

#[cfg(test)]
mod tests {
    use std::{fs, time::Instant};

    use rand::RngCore;

    use super::*;
    
    // cargo test --release --package opengm_rts --bin opengm_rts --features build-binary -- tests::test_rts --exact --show-output 
    // 1亿比特
    // other( without FFT and LinearComplexity): 125s, 12G
    // all: 716s
    #[test]
    fn test_rts() {
        // const NBITS :usize = 100*100*10000 / 8;
        const NBITS :usize = 100*10000 / 8;
        let testers = get_testers(&ALL_TESTS_FUNCS, NBITS);
        // println!("{:?}", testers);
        let mut samples: Vec<Sample> = Vec::new();
        let mut data = vec![0u8; NBITS];
        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            rng.fill_bytes(&mut data);
            samples.push(Sample::from(data.as_slice()));
        }
        let start = Instant::now();
        let _result = randomness_test(&samples, &testers);
        let elapsed = Instant::now() - start;
        println!("Used time: {} s", elapsed.as_secs_f64())
    }

    #[test]
    fn test_gen_data() {
        const NBITS :usize = 100*100*10000 / 8;
        let mut data = vec![0u8; NBITS];
        let mut rng = rand::thread_rng();
        for i in 0..1000 {
            rng.fill_bytes(&mut data);
            let filename = format!("./data1y/{:03}.bin", i);
            fs::write(filename, &data).unwrap();
        }

    }


    #[test]
    fn test_waterline() {
        assert_eq!(waterline(ALPHA, 1000), 981);
    }
}
