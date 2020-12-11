pub struct Timer {
    start_time: f64
}
pub fn get_time_sec() -> f64 {
    let t = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap();
	t.as_secs() as f64 + t.subsec_nanos() as f64 * 1e-9
}

impl Timer {
    pub fn new() -> Timer {
        Timer {start_time: get_time_sec()}
    }

    pub fn get_time(&self) -> f64 {
        get_time_sec() - self.start_time
    }
}
#[test]
fn test_timer(){
    let timer = Timer::new();

    std::thread::sleep(std::time::Duration::from_millis(500));
    let t= timer.get_time();
    println!("std: {}", t);
    assert!((t-0.5).abs() < 0.01);

    std::thread::sleep(std::time::Duration::from_millis(2000));
    let t= timer.get_time();
    println!("std: {}", t);
    assert!((t-2.5).abs() < 0.01);

    std::thread::sleep(std::time::Duration::from_millis(20));
    let t= timer.get_time();
    println!("std: {}", t);
    assert!((t-2.52).abs() < 0.01);

    std::thread::sleep(std::time::Duration::from_millis(10000));
    let t= timer.get_time();
    println!("std: {}", t);
    assert!((t-12.52).abs() < 0.01);
}