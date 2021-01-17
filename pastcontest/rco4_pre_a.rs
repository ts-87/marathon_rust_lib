#![allow(unused_imports)]
#![allow(dead_code)]

use proconio::marker;
use proconio::source::auto::AutoSource;
use proconio::input;
use std::io::{BufRead, BufReader, Read};

const TIMELIMIT: f64 = 1.97;

const N: usize = 1000;
const W: usize = 8;
const K: usize = 6;
const V: i32 = 8;
const H: usize = N / W;

struct State {
    rnd: XorShift32,
    timer: Timer,
    ln: [f64; 65536],
    bestscore: i32,
    cv: Vec<(usize,usize)>,
    grid: Vec<Vec<i32>>,
    rgrid: Vec<Vec<i32>>,
    cnt: [[i32; W]; H+1],
    rowscore: [i32; H+1]
}

fn mat_swap<T: Copy>(v: &mut Vec<Vec<T>>, i1: usize, j1: usize, i2: usize, j2: usize) {
    let t = v[i2][j2];
    v[i2][j2] = v[i1][j1];
    v[i1][j1] = t;
}
/*
fn mat_swap<T>(v: &mut Vec<Vec<T>>, i1: usize, j1: usize, i2: usize, j2: usize) {
    let p1: *mut T = &mut v[i1][j1];
    let p2: *mut T = &mut v[i2][j2];
    unsafe {
        p1.swap(p2);
    }
}
*/

fn localsearch(st: &mut State) -> Vec<Vec<i32>> {
    let mut curscore = st.bestscore;
    let mut best = st.grid.clone();
    let mut turn = 0;
    let start = st.timer.get_time();
    let starttemp = 2.5;
    let endtemp = 0.001;
    let invtl = 1.0 / (TIMELIMIT - start);
    let mut cand: Vec<usize> = Vec::new();

    loop {
        let t = st.timer.get_time();
        if t > TIMELIMIT {
            break;
        }
        let ts = starttemp + (endtemp - starttemp) * (t - start) * invtl;
        for lh in 1..=H {
            for lw in 0..W {
                for rw in 0..W {
                    if lw == rw {continue;}
                    let rhh = st.rgrid[rw].lower_bound(&st.grid[lh][lw]);
                    let mut rhl = std::cmp::min(H, rhh);
                    cand.clear();
                    while rhl > 0 && st.grid[lh][lw] < st.rgrid[rw][rhl+1] && st.grid[lh-1][lw] < st.rgrid[rw][rhl] {
                        if st.grid[lh+1][lw] > st.rgrid[rw][rhl] {
                            cand.push(rhl);
                        }
                        rhl -= 1;
                    }
                    let csz = cand.len();
                    if csz == 0 {continue;}
                    let rh = cand[st.rnd.nextn(csz)];
                    turn += 1;
                    if lh == rh {
                        mat_swap(&mut st.grid, lh, lw, rh, rw);
                        mat_swap(&mut st.rgrid, lw, lh, rw, rh);
                        continue;
                    }
                    let ifr = st.cv[st.grid[rh][rw] as usize];
                    let ifl = st.cv[st.grid[lh][lw] as usize];
                    st.cnt[lh][ifr.0] += ifr.1 as i32;
                    st.cnt[lh][ifl.0] -= ifl.1 as i32;
                    st.cnt[rh][ifl.0] += ifl.1 as i32;
                    st.cnt[rh][ifr.0] -= ifr.1 as i32;
                    
                    let mut mx1 = 0;
                    let mut mx2 = 0;
                    for i in 0..K {
                        if mx1 < st.cnt[lh][i] {mx1 = st.cnt[lh][i];}
                        if mx2 < st.cnt[rh][i] {mx2 = st.cnt[rh][i];}
                    }

                    let diff: i32 = mx1 + mx2 - st.rowscore[lh] - st.rowscore[rh];
                    if diff >=0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                        curscore += diff;
                        st.rowscore[lh] = mx1;st.rowscore[rh] = mx2;
                        mat_swap(&mut st.grid, lh, lw, rh, rw);
                        mat_swap(&mut st.rgrid, lw, lh, rw, rh);
                        if st.bestscore < curscore {
                            st.bestscore = curscore;
                            best = st.grid.clone();
                        }
                    }
                    else {
                        st.cnt[lh][ifr.0] -= ifr.1 as i32;
                        st.cnt[lh][ifl.0] += ifl.1 as i32;
                        st.cnt[rh][ifl.0] -= ifl.1 as i32;
                        st.cnt[rh][ifr.0] += ifr.1 as i32;
                    }
                }
            }
        }
    }

    eprintln!("num_iter: {}, bestscore: {}", turn, st.bestscore);
    best
}

fn main() {
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let mut iter = buffer.split_whitespace();

    //let n: usize = iter.next().unwrap().parse().unwrap();
    for _ in 0..4 {iter.next().unwrap();}
    let cv: Vec<(usize, usize)> = (0..N).map(|_| {
        (iter.next().unwrap().parse().unwrap(),
        iter.next().unwrap().parse().unwrap())
    }).collect();

    let timer = Timer::new();
    let rnd = XorShift32::new(0);

    let mut grid: Vec<Vec<i32>> = vec![vec![0; W]; H+2];
    let mut rgrid: Vec<Vec<i32>> = vec![vec![0; H+2]; W];
    let mut curh: [usize; W] = [0; W];
    let mut curw: [usize; H+1] = [0; H+1];
    let mut cntcolor = [[0; W]; H+1];
    let mut lst = 0;
    for i in 0..N {
        let mut mx: i32 = -1;
        let mut nw: usize = 0;
        let mut nh: usize = 0;
        for j in 0..W {
            if curh[j] < H && curh[j] - lst < 50 && cntcolor[curh[j]+1][cv[i].0] + cv[i].1 as i32 > mx {
                mx = cntcolor[curh[j]+1][cv[i].0] + cv[i].1 as i32;
                nw = j;
                nh = curh[j] + 1;
            }
        }
        curh[nw] += 1;
        curw[nh] += 1;
        if curw[nh] == W {lst += 1;}
        cntcolor[nh][cv[i].0] += cv[i].1 as i32;
        grid[nh][nw] = i as i32;
        rgrid[nw][nh] = i as i32;
    }
    let mut score = 0;
    let mut rowscore = [0; H+1];
    for i in 0..H {
        let mut mx = 0;
        for j in 0..K {
            if mx < cntcolor[i+1][j] {
                mx = cntcolor[i+1][j];
            }
        }
        rowscore[i+1] = mx;
        score += mx;
    }

    for i in 0..W {
        grid[0][i] = -1;
        grid[126][i] = N as i32;
        rgrid[i][0] = -1;
        rgrid[i][126] = N as i32;
    }
    
    let mut ln16: [f64; 65536] = [0.0; 65536];
    for i in 0..65536 {
        ln16[i] = (i as f64 / 65536.0 + 1.0/(2.0*65536.0)).ln();
    }
    let mut state = State{rnd: rnd, timer: timer, ln: ln16, bestscore: score, cv: cv,
                                grid: grid, rgrid: rgrid, cnt: cntcolor, rowscore: rowscore};
    
    eprintln!("greedy: {}", score);

    let bestgrid = localsearch(&mut state);

    let mut ans: [usize; 1000] = [0; 1000];
    for i in 0..N {
        ans[bestgrid[i / W +1][i % W] as usize] = i % W;
    }
    for i in 0..N {
        print!("{}\n", ans[i]);
    }
}


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

pub struct XorShift32{
    pub y: u32
}

const INV32: f64 = 1.0 / std::u32::MAX as f64;
impl XorShift32 {
    pub fn new(seed:u32) -> Self {
        Self {y: seed ^ 2463534242}
    }

    #[inline]
    pub fn next_int(&mut self) -> usize {
        self.y ^= self.y << 13;
        self.y ^= self.y >> 17;
        self.y ^= self.y << 5;
        self.y as usize
    }
    
    #[inline]
    pub fn nextn(&mut self, n: usize) -> usize {
        self.next_int() % n
    }

    #[inline]
    pub fn next_double(&mut self) -> f64 {
        self.next_int() as f64 * INV32
    }
}

pub trait BinarySearch<T> {
    fn lower_bound(&self, x: &T) -> usize;
    fn upper_bound(&self, x: &T) -> usize;
}

impl<T: Ord> BinarySearch<T> for [T] {
    fn lower_bound(&self, x: &T) -> usize {
        let mut lo = 0;
        let mut hi = self.len();

        while lo != hi {
            let mid = (lo + hi) >> 1;
            if self[mid] < *x {
                lo = mid + 1;
            }
            else {
                hi = mid;
            }
        }
        lo
    }

    fn upper_bound(&self, x: &T) -> usize {
        let mut lo = 0;
        let mut hi = self.len();

        while lo != hi {
            let mid = (lo + hi) >> 1;
            if self[mid] > *x {
                hi = mid;
            }
            else {
                lo = mid + 1;
            }
        }
        lo
    }
}