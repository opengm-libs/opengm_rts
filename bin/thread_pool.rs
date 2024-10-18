use opengm_rts::*;
use std::{
    collections::HashMap,
    sync::{mpsc, Arc, Mutex},
    thread,
};
// The ThreadPool from ch20 of the book

pub struct ThreadPool {
    #[allow(dead_code)]
    workers: Option<Vec<Worker>>,
    job_sender: mpsc::SyncSender<Job>,
    // result_receiver: mpsc::Receiver<HashMap<Tester, TestResult>>,
}

type Job = Box<dyn FnOnce() -> HashMap<Tester, TestResult> + Send + 'static>;

impl ThreadPool {
    pub fn new(size: usize) -> (ThreadPool, mpsc::Receiver<HashMap<Tester, TestResult>>) {
        assert!(size > 0);
        let num_cpus = thread::available_parallelism().unwrap().get();

        let (job_sender, job_receiver) = mpsc::sync_channel(num_cpus);
        let (result_sender, result_receiver) = mpsc::channel();

        let mut workers = Vec::with_capacity(size);

        let job_receiver = Arc::new(Mutex::new(job_receiver));
        let result_sender = Arc::new(Mutex::new(result_sender));
        for id in 0..size {
            workers.push(Worker::new(id, Arc::clone(&job_receiver), Arc::clone(&result_sender)));
        }
        (
            ThreadPool {
                workers: Some(workers),
                job_sender,
            },
            result_receiver,
        )
    }

    pub fn execute<F>(&self, f: F)
    where
        F: FnOnce() -> HashMap<Tester, TestResult> + Send + 'static,
    {
        let job = Box::new(f);
        self.job_sender.send(job).unwrap();
    }
}

// impl Drop for ThreadPool {
//     fn drop(&mut self) {
//         for _ in 0..self.workers.as_ref().unwrap().len() {
//             self.sender.send(Job::Terminate).unwrap();
//         }

//         let workers = self.workers.take().unwrap();
//         for worker in workers {
//             println!("Shutting down worker {}", worker.id);
//             worker.thread.join().unwrap();
//         }
//     }
// }

#[allow(dead_code)]
struct Worker {
    id: usize,
    thread: thread::JoinHandle<()>,
}

impl Worker {
    fn new(
        id: usize,
        job_receiver: Arc<Mutex<mpsc::Receiver<Job>>>,
        result_sender: Arc<Mutex<mpsc::Sender<HashMap<Tester, TestResult>>>>,
    ) -> Worker {
        // println!("start worker {}", id);

        let thread = thread::spawn(move || loop {
            let job = job_receiver.lock().unwrap().recv();
            match job {
                Ok(f) => {
                    let result = f();
                    result_sender.lock().unwrap().send(result).unwrap();
                }
                Err(_) => break,
            }
        });
        Worker { id, thread }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        thread,
        time::Instant,
    };

    use opengm_rts::{get_testers, sample_test, Sample, ALL_TESTS_FUNCS};
    use rand::RngCore;

    use super::ThreadPool;

    // #[test]
    // fn test_thread_pool() {
    //     const NBITS: usize = 100 * 100 * 10000 / 8;
    //     // const NBITS: usize = 100 * 10000 / 8;
    //     let testers = get_testers(&ALL_TESTS_FUNCS, NBITS);
    //     let mut samples: Vec<Sample> = Vec::new();
    //     let mut data = vec![0u8; NBITS];
    //     let mut rng = rand::thread_rng();
    //     for _ in 0..1000 {
    //         rng.fill_bytes(&mut data);
    //         samples.push(Sample::from(data.as_slice()));
    //     }

    //     let start = Instant::now();

    //     let num_cpus = thread::available_parallelism().unwrap().get();
    //     let pool = ThreadPool::new(num_cpus);

    //     for sample in samples {
    //         let sample = sample;
    //         let testers = testers.clone();
    //         pool.execute(move || sample_test(&sample, &testers));
    //     }

    //     let mut results = Vec::with_capacity(1000);
    //     for _ in 0..1000 {
    //         results.push(pool.result_receiver.recv().unwrap());
    //     }

    //     let elapsed = Instant::now() - start;
    //     println!("Used time: {} s", elapsed.as_secs_f64());

    //     println!("{}", results.len());
    // }

    #[test]
    fn test_thread_pool_small_memory() {
        // const NBITS :usize = 100*100*10000 / 8;
        const NBITS: usize = 100 * 10000 / 8;

        let testers = get_testers(&ALL_TESTS_FUNCS, NBITS);

        let start = Instant::now();

        let num_cpus = thread::available_parallelism().unwrap().get();

        // The statistics threads pool.
        let (pool, result_receiver) = ThreadPool::new(num_cpus);

        // The reading thread
        thread::spawn(move || {
            let mut rng = rand::thread_rng();
            let mut data = vec![0u8; NBITS];
            for i in 0..1000 {
                println!("sample: {}", i);
                let testers = testers.clone();
                rng.fill_bytes(&mut data);
                let sample = Sample::from(data.as_slice());
                pool.execute(move || sample_test(&sample, &testers));
            }
        });

        // In main thread, read results from pool.
        let mut results = Vec::with_capacity(1000);
        for _ in 0..1000 {
            results.push(result_receiver.recv().unwrap());
        }

        let elapsed = Instant::now() - start;
        println!("Used time: {} s", elapsed.as_secs_f64());

        println!("{}", results.len());
    }
}
