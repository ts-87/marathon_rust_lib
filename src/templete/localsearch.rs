const TIMELIMIT: f64 = 1.97;

pub struct State {
    rnd: XorShift32,
    timer: Timer,
    ln: [f64; 65536],
    bestscore: i32,
}

pub fn localsearch(st: &mut State) -> Vec<Vec<i32>> {
    let mut curscore = st.bestscore;
    let mut best = st.grid.clone();
    let mut turn = 0;
    let start = st.timer.get_time();
    let starttemp = 2.5;
    let endtemp = 0.001;
    let invtl = 1.0 / (TIMELIMIT - start);

    loop {
        let t = st.timer.get_time();
        if t > TIMELIMIT {
            break;
        }
        let ts = starttemp + (endtemp - starttemp) * (t - start) * invtl;
        

        let diff: i32;
        if diff >=0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
            if st.bestscore < curscore {
                st.bestscore = curscore;
            }
        }
    }

    eprintln!("num_iter: {}, bestscore: {}", turn, st.bestscore);
    best
}

fn main() {
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let source = AutoSource::from(buffer.as_str());

    input! {
        from source,
        _: usize,_: usize,_: usize,_: usize,
        cv: [(usize, usize); N]
    }


    let timer = Timer::new();
    let rnd = XorShift32::new(0);

    let mut score: i32 = 0;

    let mut ln16: [f64; 65536] = [0.0; 65536];
    for i in 0..65536 {
        ln16[i] = (i as f64 / 65536.0 + 1.0/(2.0*65536.0)).ln();
    }
    
    let mut state = State{rnd: rnd, timer: timer, ln: ln16, bestscore: score};
    
    eprintln!("greedy: {}", score);
}