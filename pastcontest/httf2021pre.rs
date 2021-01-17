#![allow(unused_imports)]
#![allow(dead_code)]

use std::io::Read;

const TIMELIMIT: f64 = 2.97;

const N: usize = 20;
const M: usize = 100;
const F: i32 = 15;
const DX:[i32; 4] = [-1, 0, 1, 0];
const DY:[i32; 4] = [0, -1, 0, 1];

fn move_to(ans: &mut String, px: i32, py: i32, nx: i32, ny: i32) {
        if nx > px {
            for _ in 0..(nx-px){
                ans.push('D');
            }
        }
        else{
            for _ in 0..(px-nx){
                ans.push('U');
            }
        }
        if ny > py {
            for _ in 0..(ny-py){
                ans.push('R');
            }
        }
        else{
            for _ in 0..(py-ny){
                ans.push('L');
            }
        }
}

struct State {
    rnd: XorShift32,
    timer: Timer,
    ln: [f64; 65536],
    bestscore: i32,
    pos: Vec<(i32,i32)>,
    card: Vec<usize>
}
fn tsp(st: &mut State) -> Vec<usize> {
    let mut curscore = st.bestscore;
    let mut bestscore = curscore;
    let mut best = st.card.clone();
    let mut turn = 0;
    let start = st.timer.get_time();
    let starttemp = 2.0;
    let endtemp = 0.001;
    let limit = 0.2 + start;
    let invtl = 1.0 / 0.2;
    //st.pos.push((10,10));st.card.push(M);
    let n = st.card.len();
    loop {
        let t = st.timer.get_time();
        if t > limit {
            break;
        }
        let ts = starttemp + (endtemp - starttemp) * (t - start) * invtl;
        for id1 in 1..n {
            for id2 in id1+1..n {
                turn += 1;
                /*
                let bp1 = match (st.card[id1] as i32 - st.card[id1-1] as i32).abs() {
                    1 => -F,
                    _ => 0,
                };
                let bp2 = match (st.card[id2] as i32 - st.card[id2+1] as i32).abs() {
                    1 => -F,
                    _ => 0,
                 };
                let bn1 = match (st.card[id1-1] as i32 - st.card[id2] as i32).abs() {
                    1 => -F,
                    _ => 0,
                }; 
                let bn2 = match (st.card[id1] as i32 - st.card[id2+1] as i32).abs() {
                    1 => -F,
                    _ => 0,
                };
                */
                let pre1= (st.pos[st.card[id1-1]].0 - st.pos[st.card[id1]].0).abs() + (st.pos[st.card[id1-1]].1 - st.pos[st.card[id1]].1).abs();
                let pre2= if id2 < n-1 {
                    (st.pos[st.card[id2]].0 - st.pos[st.card[id2+1]].0).abs() + (st.pos[st.card[id2]].1 - st.pos[st.card[id2+1]].1).abs()
                }
                else {0};
                let nx1= (st.pos[st.card[id1-1]].0 - st.pos[st.card[id2]].0).abs() + (st.pos[st.card[id1-1]].1 - st.pos[st.card[id2]].1).abs();
                let nx2= if id2 < n-1 {
                    (st.pos[st.card[id1]].0 - st.pos[st.card[id2+1]].0).abs() + (st.pos[st.card[id1]].1 - st.pos[st.card[id2+1]].1).abs()
                }
                else {0};

                let diff: i32 = -nx1 - nx2 + pre1 + pre2;// - bp1 - bp2 + bn1 + bn2;
                if diff >=0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                    curscore -= diff;
                    st.card[id1..id2+1].reverse();
                    if bestscore > curscore {
                        bestscore = curscore;
                        best = st.card.clone();
                    }
                }
            }
        }
    }
    st.bestscore = bestscore;
    eprintln!("num_iter: {}, bestscore: {}", turn, bestscore);
    best
}

fn sa(st: &mut State, ip: (i32, i32)) -> Vec<(i32,i32)> {
    let mut curscore = 0;
    let mut best = st.pos.clone();
    let mut turn = 0;
    let start = st.timer.get_time();
    let starttemp = 5.0;
    let endtemp = 0.001;
    let invtl = 1.0 / (TIMELIMIT - start);

    curscore += (ip.0 - st.pos[st.card[0]].0).abs() + (ip.1 - st.pos[st.card[0]].1).abs();
    for i in 1..M {
        curscore += (st.pos[st.card[i]].0 - st.pos[st.card[i-1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i-1]].1).abs();
    }
    curscore += (st.pos[st.card[M-1]].0 - st.pos[0].0).abs() + (st.pos[st.card[M-1]].1 - st.pos[0].1).abs();
    for i in 1..M {
        curscore += (st.pos[i].0 - st.pos[i-1].0).abs() + (st.pos[i].1 - st.pos[i-1].1).abs();
    }
    let mut bestscore = curscore;
    eprintln!("{}", bestscore);
    loop {
        let t = st.timer.get_time();
        if t > TIMELIMIT {
            break;
        }
        let ts = starttemp + (endtemp - starttemp) * (t - start) * invtl;
        for i in 0..M {
            for j in i+1..M {
                turn += 1;

                let mut pre = 0;
                let mut nx = 0;
                if i == 0 {
                    pre += (ip.0 - st.pos[st.card[i]].0).abs() + (ip.1 - st.pos[st.card[i]].1).abs()
                        +  (st.pos[st.card[i]].0 - st.pos[st.card[i+1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i+1]].1).abs();
                }
                else {
                    pre += (st.pos[st.card[i]].0 - st.pos[st.card[i-1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i-1]].1).abs()
                        +  (st.pos[st.card[i]].0 - st.pos[st.card[i+1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i+1]].1).abs();
                }
                if j == M-1 {
                    pre += (st.pos[st.card[j]].0 - st.pos[0].0).abs() + (st.pos[st.card[j]].1 - st.pos[0].1).abs()
                        +  (st.pos[st.card[j]].0 - st.pos[st.card[j-1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j-1]].1).abs();
                }
                else {
                    pre += (st.pos[st.card[j]].0 - st.pos[st.card[j+1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j+1]].1).abs()
                        +  (st.pos[st.card[j]].0 - st.pos[st.card[j-1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j-1]].1).abs();
                }

                if st.card[i] == 0 {
                    pre += (st.pos[st.card[i]].0 - st.pos[st.card[M-1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[M-1]].1).abs()
                        +(st.pos[st.card[i]].0 - st.pos[st.card[i]+1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]+1].1).abs();
                }
                else if st.card[i] == M-1 {
                    pre += (st.pos[st.card[i]].0 - st.pos[st.card[i]-1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]-1].1).abs();
                }
                else {
                    pre += (st.pos[st.card[i]].0 - st.pos[st.card[i]-1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]-1].1).abs()
                        +(st.pos[st.card[i]].0 - st.pos[st.card[i]+1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]+1].1).abs();
                }

                if st.card[j] == 0 {
                    pre += (st.pos[st.card[j]].0 - st.pos[st.card[M-1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[M-1]].1).abs()
                        +(st.pos[st.card[j]].0 - st.pos[st.card[j]+1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]+1].1).abs();
                }
                else if st.card[j] == M-1 {
                    pre += (st.pos[st.card[j]].0 - st.pos[st.card[j]-1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]-1].1).abs();
                }
                else {
                    pre += (st.pos[st.card[j]].0 - st.pos[st.card[j]-1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]-1].1).abs()
                        +(st.pos[st.card[j]].0 - st.pos[st.card[j]+1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]+1].1).abs();
                }
                st.pos.swap(st.card[i], st.card[j]);
                if i == 0 {
                    nx += (ip.0 - st.pos[st.card[i]].0).abs() + (ip.1 - st.pos[st.card[i]].1).abs()
                        +  (st.pos[st.card[i]].0 - st.pos[st.card[i+1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i+1]].1).abs();
                }
                else {
                    nx += (st.pos[st.card[i]].0 - st.pos[st.card[i-1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i-1]].1).abs()
                        +  (st.pos[st.card[i]].0 - st.pos[st.card[i+1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i+1]].1).abs();
                }
                if j == M-1 {
                    nx += (st.pos[st.card[j]].0 - st.pos[0].0).abs() + (st.pos[st.card[j]].1 - st.pos[0].1).abs()
                        +  (st.pos[st.card[j]].0 - st.pos[st.card[j-1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j-1]].1).abs();
                }
                else {
                    nx += (st.pos[st.card[j]].0 - st.pos[st.card[j+1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j+1]].1).abs()
                        +  (st.pos[st.card[j]].0 - st.pos[st.card[j-1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j-1]].1).abs();
                }

                if st.card[i] == 0 {
                    nx += (st.pos[st.card[i]].0 - st.pos[st.card[M-1]].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[M-1]].1).abs()
                        +(st.pos[st.card[i]].0 - st.pos[st.card[i]+1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]+1].1).abs();
                }
                else if st.card[i] == M-1 {
                    nx += (st.pos[st.card[i]].0 - st.pos[st.card[i]-1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]-1].1).abs();
                }
                else {
                    nx += (st.pos[st.card[i]].0 - st.pos[st.card[i]-1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]-1].1).abs()
                        +(st.pos[st.card[i]].0 - st.pos[st.card[i]+1].0).abs() + (st.pos[st.card[i]].1 - st.pos[st.card[i]+1].1).abs();
                }

                if st.card[j] == 0 {
                    nx += (st.pos[st.card[j]].0 - st.pos[st.card[M-1]].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[M-1]].1).abs()
                        +(st.pos[st.card[j]].0 - st.pos[st.card[j]+1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]+1].1).abs();
                }
                else if st.card[j] == M-1 {
                    nx += (st.pos[st.card[j]].0 - st.pos[st.card[j]-1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]-1].1).abs();
                }
                else {
                    nx += (st.pos[st.card[j]].0 - st.pos[st.card[j]-1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]-1].1).abs()
                        +(st.pos[st.card[j]].0 - st.pos[st.card[j]+1].0).abs() + (st.pos[st.card[j]].1 - st.pos[st.card[j]+1].1).abs();
                }
                let diff = pre - nx;
                if diff >=0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                    curscore -= diff;
                    //st.pos.swap(st.card[i], st.card[j]);
                    if bestscore > curscore {
                        bestscore = curscore;
                        best = st.pos.clone();
                    }
                }
                else {st.pos.swap(st.card[i], st.card[j]);}
            }
        }
    }

    eprintln!("num_iter: {}, bestscore: {}", turn, bestscore);
    best
}
fn main() {
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let mut iter = buffer.split_whitespace();

    let mut pos: Vec<(i32, i32)> = (0..M).map(|_| {
        (iter.next().unwrap().parse().unwrap(),
        iter.next().unwrap().parse().unwrap())
    }).collect();

    let mut grid: [[usize; N]; N] = [[M; N]; N];
    for i in 0..M {
        grid[pos[i].0 as usize][pos[i].1 as usize] = i;
    }
    let mut ans: String = String::new();

    let timer = Timer::new();
    let rnd = XorShift32::new(0);
    let mut ln16: [f64; 65536] = [0.0; 65536];
    for i in 0..65536 {
        ln16[i] = (i as f64 / 65536.0 + 1.0/(2.0*65536.0)).ln();
    }

    let mut px: i32 = 0;
    let mut py: i32 = 0;
    let mut card: Vec<usize> = Vec::new();
    let mut used = [false; M];
    //let mut pre = M+10;
    let mut score = 0;
    let mut cid =0;
    let mut cardlist: Vec<usize> = Vec::new();
    let mut state = State{rnd: rnd, timer: timer, ln: ln16, bestscore: score, pos: pos,
        card: card};
    for lp in 0..1 {
        state.card.clear();
        score = 0;
        for _ in 0..M {
            let mut mn = 100000;
            let mut nid = M;
            for i in 0..M {
                if used[i] {continue;}
                /*
                match lp {
                    0 => {
                        if i >= 33 {continue;}
                    },
                    1 => {
                        if i < 33 || i >= 64 {continue;}
                    },
                    _ => {
                        if i < 64 {continue;}
                    },
                }
                */
                let (nx, ny) = state.pos[i];
                let dist = (px - nx).abs() + (py - ny).abs();
                if mn > dist {
                    mn = dist;
                    nid = i;
                }
            }
            if nid == M {break;}
            state.card.push(nid);
            used[nid] = true;
            //move_to(&mut ans, px, py, pos[nid].0, pos[nid].1);
            px = state.pos[nid].0;
            py = state.pos[nid].1;
            grid[px as usize][py as usize] = M;
            score += mn;
            //if (nid as i32 - pre as i32).abs() == 1 {score -= F;}
            //pre = nid;
        }
        state.bestscore = score;

        let cardsz = cardlist.len();
        if cardsz == 0{
            px = 0;
            py = 0;
        }
        else {
            px = state.pos[cardlist[cardsz-1]].0;
            py = state.pos[cardlist[cardsz-1]].1;
        }
        cardlist.append(&mut tsp(&mut state).clone());

        let cardsz = cardlist.len();
        for i in cid..cardsz {
            let (nx, ny) = state.pos[cardlist[i]];
            move_to(&mut ans,px, py , nx, ny);
            ans.push('I');
            score += (px - nx).abs() + (py - ny).abs();
            px = nx;
            py = ny;
            //pos[card[i]] = (-1, -1);
        }
        cid = cardsz;
    }

    let cardsz = cardlist.len();
    cardlist.reverse();
    let rx = std::cmp::min(px + 10, 20);
    let ry = std::cmp::min(py + 10, 20);
    let lx = rx - 10;
    let ly = ry - 10;
    let mut cardid:usize = 0;
    for i in lx..rx {
        for j in ly..ry {
            state.pos[cardlist[cardid]] = (i, j);
            cardid += 1;
        }    
    }
    state.card = cardlist.clone();
    pos = sa(&mut state,(px, py));

    let mut cardid:usize = 0;
    loop {
        let (nx, ny) = pos[cardlist[cardid]];
        move_to(&mut ans, px, py, nx, ny);
        px = nx;
        py = ny;
        if grid[px as usize][py as usize] == M {
            grid[px as usize][py as usize] = cardlist[cardid];
            //pos[cardlist[cardid]] = (px, py);
            ans.push('O');
            cardid += 1;
        }
        if cardid == cardsz {break;}
    }

    for i in 0..M {
        let (nx, ny) = pos[i];
        /*
        let lx = std::cmp::min(px, nx);
        let ly = std::cmp::min(py, ny);
        let mut rx = std::cmp::max(px, nx);
        let ry = std::cmp::max(py, ny);

        let mut mnid = M;
        let mut tx = 0;
        let mut ty = 0;
        rx = std::cmp::min(14, rx);
        for j1 in ly..=ry {
            for j2 in lx..=rx {
                if grid[j2 as usize][j1 as usize] >= 50 && grid[j2 as usize][j1 as usize] < mnid && grid[j2 as usize][j1 as usize] > i {
                    mnid = grid[j2 as usize][j1 as usize];
                    ty = j1;
                    tx = j2;
                }
            }
        }

        let mut lx = std::cmp::min(tx, nx);
        let ly = std::cmp::min(ty, ny);
        let rx = std::cmp::max(tx, nx);
        let ry = std::cmp::max(ty, ny);

        lx = std::cmp::max(tx, 17);
        let mut ttx = -1;
        let mut tty = -1;
        for j1 in ly..=ry {
            for j2 in lx..=rx {
                if grid[j2 as usize][j1 as usize] == M {
                    tty = j1;
                    ttx = j2;
                }
            }
        }
        if mnid < M && ttx != -1 {
            move_to(&mut ans, px, py , tx, ty);
            ans.push('I');
            grid[tx as usize][ty as usize] = M;
            move_to(&mut ans,tx, ty , ttx, tty);
            ans.push('O');
            grid[ttx as usize][tty as usize] = mnid;
            px = ttx;
            py = tty;
            pos[mnid]=(ttx, tty);
            
        }
        */
        move_to(&mut ans, px, py , nx, ny);
        ans.push('I');
        px = nx;
        py = ny;
        grid[px as usize][py as usize] = M;
    }
    
    println!("{}", ans);
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

