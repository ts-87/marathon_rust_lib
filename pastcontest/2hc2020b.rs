#![allow(unused_imports)]
#![allow(dead_code)]

use std::{collections::{VecDeque, HashSet}, io::{StdoutLock, prelude::*}, writeln};
use std::io::{BufWriter};

const OUT_INFO: bool = false;
//interactive io macro
//https://atcoder.jp/contests/hokudai-hitachi2019-1/submissions/8550933
pub fn readln() -> String {
	let mut line = String::new();
	::std::io::stdin().read_line(&mut line).unwrap_or_else(|e| panic!("{}", e));
	line
}
 
macro_rules! read {
	($($t:tt),*; $n:expr) => {{
		let stdin = ::std::io::stdin();
		let ret = ::std::io::BufRead::lines(stdin.lock()).take($n).map(|line| {
			let line = line.unwrap();
			let mut it = line.split_whitespace();
			_read!(it; $($t),*)
		}).collect::<Vec<_>>();
		ret
	}};
	($($t:tt),*) => {{
		let line = readln();
		let mut it = line.split_whitespace();
		_read!(it; $($t),*)
	}};
}
 
macro_rules! _read {
	($it:ident; [char]) => {
		_read!($it; String).chars().collect::<Vec<_>>()
	};
	($it:ident; [u8]) => {
		Vec::from(_read!($it; String).into_bytes())
	};
	($it:ident; usize1) => {
        $it.next().unwrap_or_else(|| panic!("input mismatch")).parse::<usize>().unwrap_or_else(|e| panic!("{}", e)) - 1
	};
	($it:ident; [usize1]) => {
		$it.map(|s| s.parse::<usize>().unwrap_or_else(|e| panic!("{}", e)) - 1).collect::<Vec<_>>()
	};
	($it:ident; [$t:ty]) => {
		$it.map(|s| s.parse::<$t>().unwrap_or_else(|e| panic!("{}", e))).collect::<Vec<_>>()
	};
	($it:ident; $t:ty) => {
		$it.next().unwrap_or_else(|| panic!("input mismatch")).parse::<$t>().unwrap_or_else(|e| panic!("{}", e))
	};
	($it:ident; $($t:tt),+) => {
		($(_read!($it; $t)),*)
	};
}

const TIMELIMIT: f64 = 60.;

const BASESCORE: f64 = 3000000.;
const V: usize = 225;
const T_MAX: usize = 1000;
const GAMMA: f64 = 2.0;
const DIV_RANGE: usize = 50;
const INF: usize = std::usize::MAX / 3;
const INF32: u32 = std::u32::MAX / 3;
const V8: u8 = 225;
const CHARGE_FROM: usize = 1<<15;
const CHARGE_TO: usize = 1<<16;

enum DayType {
    Fine,
    FineSquall,
    Rainy,
    RainyFine
}

struct Graph {
    v: usize,
    e: usize,
    graph: Vec<Vec<usize>>,
    costs: Vec<Vec<usize>>,
    next: Vec<Vec<usize>>,
    nearest_grid: Vec<(usize, usize)>,
}

struct NanoGrids {
    n_grid: usize, n_div: usize,
    p_event: f64, delta_event: i32,
    n_pattern: usize, sigma_ele: usize,
    pw_predicts: Vec<Vec<i32>>,
    c_grid_init: usize, c_grid_max: usize, v_grid_max: usize,
    x: Vec<usize>,
    patterns: Vec<usize>,
    day_type: DayType,

    charge: Vec<usize>,
    pos_to_id:Vec<usize>,
    pw_expected: Vec<Vec<i32>>,
    visit_ev_id: Vec<VecDeque<usize>>,
    chrg_end_time: Vec<VecDeque<usize>>,
    use_count: Vec<usize>,
    in_sudden: Vec<bool>,
    overflow: Vec<i32>,
    need: Vec<i32>,
    fsthalf_minexpect_divt: Vec<usize>,
    maxexpect_divt: Vec<usize>
}

struct EVs {
    n_ev: usize,
    c_ev_init: usize, c_ev_max: usize, v_ev_max: usize,
    n_trans_max: usize,
    delta_ev_move: usize,
    charge: Vec<usize>,
    pos: Vec<usize>,
    next: Vec<(usize, usize)>,//id,dist
    adj: Vec<Vec<usize>>,
    commands: Vec<VecDeque<Command>>,
    pre_commands: Vec<Command>,
    order: Vec<Vec<usize>>,
    order_cnt: Vec<usize>,
    to_charge_flag: Vec<bool>,
}

struct Orders {
    n_order: usize,
    ids: Vec<usize>,
    term_pos: Vec<(usize, usize)>,
    state: Vec<usize>,
    time: Vec<usize>,
    assign_ev: Vec<usize>,
}

#[derive(Clone)]
enum Command {
    Stay,
    Move(usize),
    ChargeFromGrid(usize),
    ChargeToGrid(usize),
    Pickup(usize),
    Check,
}

struct State {
    rnd: XorShift32,
    timer: Timer,
    ln: [f64; 65536],
    n_solution: usize,
    graph: Graph,
    nano_grid: NanoGrids,
    ev: EVs,
    order: Orders,
    cur_solution_id: usize,
    sudden: bool,
}

fn init_state(st: &mut State) {
    for i in st.nano_grid.charge.iter_mut() {*i = st.nano_grid.c_grid_init;}
    for gid in 0..st.nano_grid.n_grid {
        st.nano_grid.visit_ev_id[gid].clear();
        st.nano_grid.chrg_end_time[gid].clear();
        st.nano_grid.need[gid] = 0;
        st.nano_grid.in_sudden[gid] = false;
    }
    for eid in 0..st.ev.n_ev {
        st.ev.pre_commands[eid] = Command::Stay.clone();
        st.ev.commands[eid].clear();
        st.ev.order_cnt[eid] = 0;
        st.ev.to_charge_flag[eid] = false;
    }
    for i in st.order.assign_ev.iter_mut() {*i = V;}
    st.ev.pre_commands = vec![Command::Stay; st.ev.n_ev];
    for i in 0..st.ev.n_ev {
        st.ev.commands[i].clear();
        st.ev.charge[i] = st.ev.c_ev_init;
        st.ev.order[i].clear();
    }

    st.order.ids.clear();
    for i in 0..1024 {
        st.order.term_pos[i] = (0, 0);
        st.order.state[i] = 0;
        st.order.time[i] = 0;
    }
}

fn print_status(st: &State, t: usize) {
    eprint!("Turn {}\n", t);
    eprint!("EV charge  : ");
    for &i in st.ev.charge.iter() {
        eprint!("{:>5} ", i);
    }
    eprint!("\nGrid charge: ");
    for i in st.nano_grid.charge.iter() {
        eprint!("{:>5} ", i);
    }
    eprint!("\nOrder id/time: ");
    for &i in st.order.ids.iter() {
        eprint!("[{:>4}/{:>4}] ", i, st.order.time[i]);
    }
    eprint!("\n");
}

#[inline]
fn output_command(com: &Command, out: &mut BufWriter<StdoutLock>) {
    match com {
        Command::Stay => {write!(out, "stay\n").unwrap();},
        Command::Move(v) => {write!(out, "move {}\n", v+1).unwrap();},
        Command::ChargeFromGrid(c) => {write!(out, "charge_from_grid {}\n", c).unwrap();},
        Command::ChargeToGrid(c) => {write!(out, "charge_to_grid {}\n", c).unwrap();},
        Command::Pickup(w) => {write!(out, "pickup {}\n", w+1).unwrap();},
        _ => {}
    }
}

fn read_info(st: &mut State, div_t: usize) {
    let is_sudden_daytype = match st.nano_grid.day_type {
        DayType::FineSquall | DayType::RainyFine => {true},
        _ => {false}
    };
    st.sudden = false;
    for i in 0..st.nano_grid.n_grid {
        let (xi, yi, pw_actual, _pw_exceed, _pw_buy) = read!(usize1, usize, i32, usize, usize);
        let pat = st.nano_grid.patterns[i];
        if let DayType::RainyFine = st.nano_grid.day_type {st.nano_grid.in_sudden[i] = false;}
        if is_sudden_daytype && (st.nano_grid.pw_predicts[pat][div_t] - pw_actual).abs() > 600 {
            st.sudden = true;
            st.nano_grid.in_sudden[i] = true;
        }
        let id = st.nano_grid.pos_to_id[xi];
        st.nano_grid.charge[id] = yi;
    }

    for i in 0..st.ev.n_ev {
        st.ev.charge[i] = read!(usize);
        let (ui, vi, dist_ui, dist_vi) = read!(usize1, usize1, usize, usize);
        if ui == vi {
            st.ev.pos[i] = ui;
            st.ev.next[i].0 = ui;
            st.ev.next[i].1 = 0;
        }
        else {
            if st.ev.next[i].1 == 0 {
                if st.ev.next[i].0 == ui {
                    st.ev.next[i].0 = vi;
                    st.ev.next[i].1 = dist_vi;
                    
                }
                else {
                    st.ev.next[i].0 = ui;
                    st.ev.next[i].1 = dist_ui;
                }
            }
            else {
                if st.ev.next[i].0 == vi {st.ev.next[i].1 = dist_vi;}
                else {st.ev.next[i].1 = dist_ui;}
            }
        }
        let mut buffer = String::new();
        std::io::stdin().read_line(&mut buffer).unwrap();
        let mut input = buffer.split_whitespace();
        let n_adj = input.next().unwrap().parse().unwrap();
        st.ev.adj[i].resize(n_adj, 0);
        for j in 0..n_adj {
            st.ev.adj[i][j] = input.next().unwrap().parse::<usize>().unwrap() - 1;
        }

        let mut buffer = String::new();
        std::io::stdin().read_line(&mut buffer).unwrap();
        let mut input = buffer.split_whitespace();
        let n_order = input.next().unwrap().parse().unwrap();
        st.ev.order[i].resize(n_order, 0);
        for j in 0..n_order {
            st.ev.order[i][j] = input.next().unwrap().parse::<usize>().unwrap() - 1;
        }
    }
    let n_order = read!(usize);
    st.order.n_order = n_order;
    st.order.ids.clear();
    for _ in 0..n_order {
        let (id, wi, zi, sti, ti) = read!(usize1, usize1, usize1, usize, usize);
        st.order.ids.push(id);
        st.order.term_pos[id] = (wi, zi);
        st.order.state[id] = sti;
        st.order.time[id] = ti;
    }
}

    /*
    fine
    -207,-169,-123,-17,45,159,186,240,248,244,230,171,109,62,-55,-99,-133,-125,-57,-28
    48,83,128,169,206,133,122,89,96,93,121,100,72,38,-2,-126,-169,-157,-134,-110
    -116,-76,-75,-49,-12,14,39,64,88,100,90,60,-3,-38,-75,-74,-117,-100,-83,-78
    25000 14650 6200 50 -800 1450 9400 18700 30700 43100 55300 66800 75350 80800 83900 81150 76200 69550 63300 60450 59050 
    25000 27400 31550 37950 46400 56700 63350 69450 73900 78700 83350 89400 94400 98000 99900 99800 93500 85050 77200 70500 65000 
    25000 19200 15400 11650 9200 8600 9300 11250 14450 18850 23850 28350 31350 31200 29300 25550 21850 16000 11000 6850 2950 

    rainy
    67,76,87,80,113,122,91,92,96,94,50,47,8,-12,-97,-106,-105,-49,-21,8
    2,9,18,26,33,56,27,-18,-17,-17,-20,-24,11,4,-4,-13,-21,-1,22,46
    -82,-70,-65,-72,-24,-23,-16,-4,-24,-10,70,56,15,8,3,-1,-10,16,-7,-2
    25000 28350 32150 36500 40500 46150 52250 56800 61400 66200 70900 73400 75750 76150 75550 70700 65400 60150 57700 56650 57050 
    25000 25100 25550 26450 27750 29400 32200 33550 32650 31800 30950 29950 28750 29300 29500 29300 28650 27600 27550 28650 30950 
    25000 20900 17400 14150 10550 9350 8200 7400 7200 6000 5500 9000 11800 12550 12950 13100 13050 12550 13350 13000 12900 
    */

fn calc_expected_charge(st: &mut State, available_charge: &mut Vec<usize>, turn: usize) {
    let div_t = turn / DIV_RANGE;
    let mod_t = (turn % DIV_RANGE) as i32;

    match st.cur_solution_id {
        0..=2 => {
            for i in 0..st.nano_grid.n_grid {
                let pat = st.nano_grid.patterns[i];
                let mut tmp = st.nano_grid.charge[i] as i32;
                if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
                    tmp += st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]]
                    - st.nano_grid.pw_expected[pat][div_t];
                }
                else if div_t < st.nano_grid.maxexpect_divt[pat] {
                    tmp += st.nano_grid.pw_expected[pat][st.nano_grid.maxexpect_divt[pat]]
                    - st.nano_grid.pw_expected[pat][div_t];
                }
                else {
                    tmp += st.nano_grid.pw_expected[pat][st.nano_grid.n_div]
                    - st.nano_grid.pw_expected[pat][div_t];
                }

                tmp -= st.nano_grid.pw_predicts[pat][div_t] * mod_t;
                let sub = match st.nano_grid.day_type {
                    DayType::Fine => {st.nano_grid.c_grid_max},
                    DayType::Rainy => {45000},
                    DayType::FineSquall => {45000},
                    DayType::RainyFine => {st.nano_grid.c_grid_max}
                };
                //available_charge[i] = std::cmp::max(0, tmp - (st.nano_grid.c_grid_max as i32)) as usize;
                available_charge[i] = std::cmp::max(0, tmp - (sub as i32)) as usize;
                available_charge[i] = available_charge[i].min(st.nano_grid.charge[i]);
            }
        },
        3..=4 => {
            for i in 0..st.nano_grid.n_grid {
                let pat = st.nano_grid.patterns[i];
                let mut tmp = st.nano_grid.charge[i] as i32;
                if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
                    tmp += st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]]
                    - st.nano_grid.pw_expected[pat][div_t];
                }
                else {
                    tmp += st.nano_grid.pw_expected[pat][st.nano_grid.n_div]
                    - st.nano_grid.pw_expected[pat][div_t];
                    
                }
                tmp -= st.nano_grid.pw_predicts[pat][div_t] * mod_t;
                let sub = match st.nano_grid.day_type {
                    DayType::Fine => {10000},
                    DayType::Rainy => {250},
                    DayType::FineSquall => {250},
                    DayType::RainyFine => {
                        if st.cur_solution_id == 4 && st.ev.n_ev <= 16 {5000} else {10000}
                    }
                };
                //available_charge[i] = std::cmp::max(250, tmp) as usize - 250;
                available_charge[i] = std::cmp::max(0, tmp - (sub as i32)) as usize;
                available_charge[i] = available_charge[i].min(st.nano_grid.charge[i]);
            }
        },
        _ => {}
    }
    
    match st.nano_grid.day_type {
        DayType::FineSquall => {
            for gid in 0..st.nano_grid.n_grid {
                available_charge[gid] -= 5000.min(available_charge[gid]);
            }
        },
        _ => {}
    }
    
}

fn calc_need_charge(st: &mut State, turn: usize) {
    let div_t = turn / DIV_RANGE;
    let mod_t = turn % DIV_RANGE;
    for gid in 0..st.nano_grid.n_grid {
        let sudden_delta = if st.nano_grid.in_sudden[gid] {
            match st.nano_grid.day_type {
                DayType::FineSquall => {-1000},
                DayType::RainyFine => {1000},
                _ => {0}
            }
        }
        else {0};
        let pat = st.nano_grid.patterns[gid];
        let mut tmp = st.nano_grid.charge[gid] as i32;
        if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
            tmp += st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]]
            - st.nano_grid.pw_expected[pat][div_t];
        }
        else {
            tmp += st.nano_grid.pw_expected[pat][st.nano_grid.n_div]
            - st.nano_grid.pw_expected[pat][div_t];
        }
        tmp -= st.nano_grid.pw_predicts[pat][div_t] * mod_t as i32;
        //tmp += st.nano_grid.pw_expected[pat][20.min(div_t+2)] - st.nano_grid.pw_expected[pat][div_t];
        tmp += sudden_delta * ((DIV_RANGE - mod_t) as i32);
        st.nano_grid.need[gid] = std::cmp::max(0, -tmp+250);

        let mut tmp = st.nano_grid.charge[gid] as i32;
        
        if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
            tmp += st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]]
            - st.nano_grid.pw_expected[pat][div_t];
        }
        else if div_t < st.nano_grid.maxexpect_divt[pat] {
            tmp += st.nano_grid.pw_expected[pat][st.nano_grid.maxexpect_divt[pat]]
            - st.nano_grid.pw_expected[pat][div_t];
        }
        else {
            tmp += st.nano_grid.pw_expected[pat][st.nano_grid.n_div]
            - st.nano_grid.pw_expected[pat][div_t];
        }
        tmp -= st.nano_grid.pw_predicts[pat][div_t] * mod_t as i32;
        tmp += sudden_delta * ((DIV_RANGE - mod_t) as i32);
        st.nano_grid.overflow[gid] = std::cmp::max(0, tmp - st.nano_grid.c_grid_max as i32);
        st.nano_grid.overflow[gid] = std::cmp::min(st.nano_grid.charge[gid] as i32, st.nano_grid.overflow[gid]);
    }
}

fn optimaize_init_matching(st: &mut State, available_charge: &Vec<usize>) -> Vec<usize> {
    let mut charge_ev_cnt = vec![0; st.nano_grid.n_grid];
    let charge_amount_per_ev = 5000;//available grid数とev数で決める
    for i in 0..st.nano_grid.n_grid {
        charge_ev_cnt[i] = (available_charge[i] + 1000) / charge_amount_per_ev;
        charge_ev_cnt[i] = charge_ev_cnt[i].min(2);
    }
    let mut curscore = 0;
    let mut match_ev_grid = vec![V; st.ev.n_ev];
    for ct in (0..2).rev() {
        for i in 0..st.ev.n_ev {
            if match_ev_grid[i] != V {continue;}
            let mut mn = INF;
            let mut nid = V;
            for gid in 0..st.nano_grid.n_grid {
                if charge_ev_cnt[gid] <= ct {continue;}
                if mn > st.graph.costs[st.ev.next[i].0][st.nano_grid.x[gid]] {
                    mn = st.graph.costs[st.ev.next[i].0][st.nano_grid.x[gid]];
                    nid = gid;
                }
            }
            if nid != V {
                match_ev_grid[i] = st.nano_grid.x[nid];
                charge_ev_cnt[nid] -= 1;
                curscore += mn;
            }
        }
    }
    if OUT_INFO {eprint!("init_match: {}\n", curscore);}

    let mut bestscore = curscore;
    let mut best = match_ev_grid.clone();
    //let start = st.timer.get_time();
    let starttemp = 1.5;
    let endtemp = 0.001;
    let num_itr = 10000;
    let invtl = 1.0 / num_itr as f64;
    for lp in 0..num_itr {
        let ts = starttemp + (endtemp - starttemp) * lp as f64 * invtl;

        for ev1 in 0..st.ev.n_ev {
            for ev2 in 0..st.ev.n_ev {
                if ev1 == ev2 {continue;}
                let pre = st.graph.costs[st.ev.next[ev1].0][match_ev_grid[ev1]] + st.graph.costs[st.ev.next[ev2].0][match_ev_grid[ev2]];
                let nx =  st.graph.costs[st.ev.next[ev1].0][match_ev_grid[ev2]] + st.graph.costs[st.ev.next[ev2].0][match_ev_grid[ev1]];
                let diff = pre as i32 - nx as i32;
                if diff >= 0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                    curscore -= diff as usize;
                    match_ev_grid.swap(ev1, ev2);
                    if bestscore > curscore {
                        bestscore = curscore;
                        best = match_ev_grid.clone();
                    }
                }
            }
        }
    }

    if OUT_INFO {
        let mut score = 0;
        for eid in 0..st.ev.n_ev {
            score += st.graph.costs[st.ev.next[eid].0][best[eid]];
        }
        eprint!("calc: {}/ best_match: {}\n", score, bestscore);
    }
    best
}

fn optimize_ev_route(st: &mut State, perm: &mut Vec<Vec<(usize, usize)>>, num_add_order: &Vec<(usize, usize)>) -> Vec<Vec<(usize, usize)>> {
    
    let mut best = vec![Vec::new(); st.ev.n_ev];
    let mut evcost = vec![0; st.ev.n_ev];
    //init tsp dp
    let mut dp = vec![vec![INF32; 8]; 256];
    let mut pre = vec![vec![V8; 8]; 256];
    if OUT_INFO {
        let mut score = 0;
        for ev in 0..perm.len() {
            for i in 0..perm[ev].len()-1 {
                score += st.graph.costs[perm[ev][i].0][perm[ev][i+1].0];
            }
        }
        eprintln!("init: {}", score);
    }
    let mut bestscore = 0;
    for ev in 0..perm.len() {
        let m = perm[ev].len() - 1;
        if m == 0 {
            best[ev].push(perm[ev][0]);
            best[ev].push((V, INF));
            continue;
        }
        let n = 1 << m;
        let b = num_add_order[ev].0;
        let init_pos = perm[ev][0].0;
        for i in 0..m {
            if i+1 < b || ((i+1-b) & 1 == 0) {dp[1<<i][i] = st.graph.costs[init_pos][perm[ev][i+1].0] as u32;}
        }
        for k in 1..n {
            for i in 0..m {
                if ((k>>i) & 1 == 0) || dp[k][i] == INF32 {continue;}
                for j in 0..m {
                    if (k>>j) & 1 == 1 {continue;}
                    if j+1 >= b && (j+1-b) & 1 == 1 && (k>>(j-1)) & 1 == 0 {continue;}
                    if dp[k|(1<<j)][j] > dp[k][i] + st.graph.costs[perm[ev][i+1].0][perm[ev][j+1].0] as u32 {
                        dp[k|(1<<j)][j] = dp[k][i] + st.graph.costs[perm[ev][i+1].0][perm[ev][j+1].0] as u32;
                        pre[k|(1<<j)][j] = i as u8;
                    }
                }
            }
        }

        let mut cur = V;
        let mut curstate = (1<<m)-1;
        let mut mn = INF32;
        for i in 0..m {
            if dp[curstate][i] < mn {
                mn = dp[curstate][i];
                cur = i;
            } 
        }
        evcost[ev] = mn as usize;
        bestscore += mn as usize;

        best[ev].push(perm[ev][cur+1]);
        curstate ^= 1<<cur;
        cur = pre[curstate|(1<<cur)][cur] as usize;
        while cur != V {
            best[ev].push(perm[ev][cur+1]);
            curstate ^= 1<<cur;
            cur = pre[curstate|(1<<cur)][cur] as usize;
        }
        best[ev].push(perm[ev][0]);
        best[ev].reverse();
        best[ev].push((V, INF));
        for i in 0..256 {
            for j in 0..8 {
                dp[i][j] = INF32;
                pre[i][j] = V8;
            }
        }  
    }
    //ls
    if OUT_INFO {
        let mut score = 0;
        for ev in 0..best.len() {
            for i in 0..best[ev].len()-1 {
                score += st.graph.costs[best[ev][i].0][best[ev][i+1].0];
            }
        }
        eprintln!("before_ls: {} / {}", score, bestscore);
    }
    
    let mut used_order = vec![false; 1024];
    for eid in 0..st.ev.n_ev {
        for &(_, oid) in perm[eid].iter().skip(1) {
            used_order[oid] = true;
        }
    }
    let mut rem_order = Vec::new();
    for &oid in st.order.ids.iter() {
        if !used_order[oid] {rem_order.push(oid);}
    }
    let mut perm = best.clone();
    let mut curscore = bestscore;
    let starttemp = 5.0;//-7.5
    let endtemp = 0.001;
    
    let num_itr = 45000;

    let invtl = 1.0 / num_itr as f64;
    let mut tmp_set = HashSet::new();
    for lp in 0..num_itr {
        let ts = starttemp + (endtemp - starttemp) * lp as f64 * invtl;
        for eid1 in 0 ..st.ev.n_ev {
            if perm[eid1].len() <= 3 {continue;}
            let id1 = st.rnd.nextn(perm[eid1].len() - 2) + 1;
            let l_flag = if st.order.state[perm[eid1][id1].1] == 1 {true} else {false};
            let sel = st.rnd.next_int() & 1;
            if sel == 0 || l_flag {//perm内swap
                let prev = evcost[eid1];
                let mut id2 = st.rnd.nextn(perm[eid1].len() - 2) + 1;
                if id2 == id1 {
                    id2 = id2 % (perm[eid1].len() - 2) + 1;
                }
                perm[eid1].swap(id1, id2);
                let mut nx = 0;
                tmp_set.clear();
                let mut pre = V;
                for v in perm[eid1].iter().rev().skip(1) {
                    if v.1 == INF {
                        nx += st.graph.costs[pre][v.0];
                        continue;
                    }
                    if tmp_set.contains(&v.1) {
                        nx += st.graph.costs[pre][st.order.term_pos[v.1].0];
                        pre = st.order.term_pos[v.1].0;
                    } 
                    else {
                        nx += st.graph.costs[pre][st.order.term_pos[v.1].1];
                        tmp_set.insert(v.1);
                        pre = st.order.term_pos[v.1].1;
                    }
                }
                let diff = prev as i32 - nx as i32;
                if diff >= 0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                    curscore = ((curscore as i32) - diff) as usize;
                    evcost[eid1] = nx;
                    tmp_set.clear();
                    for v in perm[eid1].iter_mut().rev().skip(1) {
                        if v.1 == INF {continue;}
                        if tmp_set.contains(&v.1) {
                            *v = (st.order.term_pos[v.1].0, v.1);
                        } 
                        else {
                            *v = (st.order.term_pos[v.1].1, v.1);
                            tmp_set.insert(v.1);
                        }
                    }
                    if bestscore > curscore {
                        bestscore = curscore;
                        best = perm.clone();
                    }
                }
                else {
                    perm[eid1].swap(id1, id2);
                }
            }
            else {//perm間swap
                let mut eid2 = st.rnd.nextn(st.ev.n_ev);
                let mut ok = false;
                for _ in 0..10 {
                    if eid2 != eid1 && perm[eid2].len() > 3 {
                        ok = true;
                        break;
                    }
                    eid2 = st.rnd.nextn(st.ev.n_ev);
                }
                if !ok {continue;}
                ok = false;
                let mut id2 = st.rnd.nextn(perm[eid2].len() - 2) + 1;
                for _ in 0..10 {
                    if st.order.state[perm[eid2][id2].1] == 0 {
                        ok = true;
                        break;
                    }
                    id2 = st.rnd.nextn(perm[eid2].len() - 2) + 1;
                }
                if !ok {continue;}
                
                let oid1 = perm[eid1][id1].1;
                let pos1 = st.order.term_pos[oid1];
                let oid2 = perm[eid2][id2].1;
                let pos2 = st.order.term_pos[oid2];
                let prev = evcost[eid1] + evcost[eid2];

                let mut ins_id1 = (V, V);
                let mut ins_id2 = (V, V);
                let (mut mn1, mut mn2) = (INF, INF);
                let mut tmp_cost;
                let mut prei = perm[eid1][0].0;
                let mut prej = prei;
                
                for i in 1..perm[eid1].len() {
                    if perm[eid1][i].1 == oid1 {continue;}
                    for j in i..perm[eid1].len() {
                        if perm[eid1][j].1 == oid1 {continue;}
                        if i == j {
                            tmp_cost = st.graph.costs[prei][pos2.0] + st.graph.costs[pos2.0][pos2.1] + st.graph.costs[perm[eid1][j].0][pos2.1];
                        }
                        else {
                            tmp_cost = 
                            st.graph.costs[prei][pos2.0] + st.graph.costs[perm[eid1][i].0][pos2.0]
                            + st.graph.costs[prej][pos2.1] + st.graph.costs[perm[eid1][j].0][pos2.1];
                        }
                        prej = perm[eid1][j].0;
                        if mn1 > tmp_cost {
                            mn1 = tmp_cost;
                            ins_id1 = (i, j);
                        }
                    }
                    prei =  perm[eid1][i].0;
                    
                }
                
                let mut prei = perm[eid2][0].0;
                let mut prej = prei;
                for i in 1..perm[eid2].len() {
                    if perm[eid2][i].1 == oid2 {continue;}
                    for j in i..perm[eid2].len() {
                        if perm[eid2][j].1 == oid2 {continue;}
                        if i == j {
                            tmp_cost = st.graph.costs[prei][pos1.0] + st.graph.costs[pos1.0][pos1.1] + st.graph.costs[perm[eid2][j].0][pos1.1];
                        }
                        else {
                            tmp_cost = 
                            st.graph.costs[prei][pos1.0] + st.graph.costs[perm[eid2][i].0][pos1.0]
                            + st.graph.costs[prej][pos1.1] + st.graph.costs[perm[eid2][j].0][pos1.1];
                        }
                        prej = perm[eid2][j].0;
                        if mn2 > tmp_cost {
                            mn2 = tmp_cost;
                            ins_id2 = (i, j);
                        }
                    }
                    prei =  perm[eid2][i].0;
                }

                let mut nx1 = 0;
                let mut pre = perm[eid1][0].0;
                for i in 1..perm[eid1].len() {
                    if ins_id1.0 == i {
                        nx1 += st.graph.costs[pre][st.order.term_pos[oid2].0];
                        pre = pos2.0;
                    }
                    if ins_id1.1 == i {
                        nx1 += st.graph.costs[pre][st.order.term_pos[oid2].1];
                        pre = pos2.1;
                    }
                    if perm[eid1][i].1 != oid1 {
                        nx1 += st.graph.costs[pre][perm[eid1][i].0];
                        pre = perm[eid1][i].0;
                    }
                }
                let mut nx2 = 0;
                pre = perm[eid2][0].0;
                for i in 1..perm[eid2].len() {
                    if ins_id2.0 == i {
                        nx2 += st.graph.costs[pre][st.order.term_pos[oid1].0];
                        pre = pos1.0;
                    }
                    if ins_id2.1 == i {
                        nx2 += st.graph.costs[pre][st.order.term_pos[oid1].1];
                        pre = pos1.1;
                    }
                    if perm[eid2][i].1 != oid2 {
                        nx2 += st.graph.costs[pre][perm[eid2][i].0];
                        pre = perm[eid2][i].0;
                    }
                }
                let diff = prev as i32 - nx1 as i32 - nx2 as i32;
                if diff >= 0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                    curscore = ((curscore as i32) - diff) as usize;

                    perm[eid1].insert(ins_id1.1,(pos2.1, oid2));
                    perm[eid1].insert(ins_id1.0,(pos2.0, oid2));
                    perm[eid2].insert(ins_id2.1,(pos1.1, oid1));
                    perm[eid2].insert(ins_id2.0,(pos1.0, oid1));
                    perm[eid1] = perm[eid1].iter().filter(|x| x.1 != oid1).cloned().collect();
                    perm[eid2] = perm[eid2].iter().filter(|x| x.1 != oid2).cloned().collect();

                    evcost[eid1] = nx1;
                    evcost[eid2] = nx2;
                    if bestscore > curscore {
                        bestscore = curscore;
                        best = perm.clone();
                    }
                }
            }

            if rem_order.len() == 0 {continue;}
                       
            for _ in 0..2 {
                
                let mut id1 = st.rnd.nextn(perm[eid1].len() - 2) + 1;
                let mut ok = false;
                for _ in 0..10 {
                    if st.order.state[perm[eid1][id1].1] == 0 {
                        ok = true;
                        break;
                    }
                    id1 = st.rnd.nextn(perm[eid1].len() - 2) + 1;
                }
                if !ok {break;}
                let oid1 = perm[eid1][id1].1;

                let id2 = st.rnd.nextn(rem_order.len());
                let oid2 = rem_order[id2];
                let pos2 = st.order.term_pos[oid2];

                let mut ins_id1 = (V, V);
                let mut mn1 = INF;
                let mut tmp_cost;
                let mut prei = perm[eid1][0].0;
                let mut prej = prei;
                for i in 1..perm[eid1].len() {
                    if perm[eid1][i].1 == oid1 {continue;}
                    for j in i..perm[eid1].len() {
                        if perm[eid1][j].1 == oid1 {continue;}
                        if i == j {
                            tmp_cost = st.graph.costs[prei][pos2.0] + st.graph.costs[pos2.0][pos2.1] + st.graph.costs[perm[eid1][j].0][pos2.1];
                        }
                        else {
                            tmp_cost = 
                            st.graph.costs[prei][pos2.0] + st.graph.costs[perm[eid1][i].0][pos2.0]
                            + st.graph.costs[prej][pos2.1] + st.graph.costs[perm[eid1][j].0][pos2.1];
                        }
                        prej =  perm[eid1][j].0;
                        if mn1 > tmp_cost {
                            mn1 = tmp_cost;
                            ins_id1 = (i, j);
                        }
                    }
                    prei =  perm[eid1][i].0;
                }
                let mut nx1 = 0;
                let mut pre = perm[eid1][0].0;
                for i in 1..perm[eid1].len() {
                    if ins_id1.0 == i {
                        nx1 += st.graph.costs[pre][st.order.term_pos[oid2].0];
                        pre = pos2.0;
                    }
                    if ins_id1.1 == i {
                        nx1 += st.graph.costs[pre][st.order.term_pos[oid2].1];
                        pre = pos2.1;
                    }
                    if perm[eid1][i].1 != oid1 {
                        nx1 += st.graph.costs[pre][perm[eid1][i].0];
                        pre = perm[eid1][i].0;
                    }
                }
                let prev = evcost[eid1] as i32;
                let diff = prev as i32 - nx1 as i32;
                if diff >= 0 || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                    curscore = ((curscore as i32) - diff) as usize;
                    evcost[eid1] = nx1;
                    perm[eid1].insert(ins_id1.1,(pos2.1, oid2));
                    perm[eid1].insert(ins_id1.0,(pos2.0, oid2));
                    perm[eid1] = perm[eid1].iter().filter(|&x| x.1 != oid1).cloned().collect();
                    rem_order[id2] = oid1;
                    if bestscore > curscore {
                        bestscore = curscore;
                        best = perm.clone();
                    }
                }
            }
        } 
    }
    
    if OUT_INFO {
        let mut score = 0;
        for ev in 0..best.len() {
            for i in 0..best[ev].len()-1 {
                score += st.graph.costs[best[ev][i].0][best[ev][i+1].0];
            }
        }
        eprintln!("after_ls: {} / {}", score, bestscore);
        /*
        for ev in 0..best.len() {
            eprint!("{} ", best[ev].len());
        }
        eprint!("\n");
        for ev in 0..best.len() {
            for i in 1..best[ev].len() {
                let flag = if best[ev][i].0 == st.order.term_pos[best[ev][i].1].0 {0} else {1};
                eprint!("{}/{} ", best[ev][i].1,flag);
            }
            eprint!("\n");
        }
        */
    }
    
    best
}

fn construct_ev_route(st: &mut State, perm: &mut Vec<Vec<(usize, usize)>>, _turn: usize) -> Vec<VecDeque<(usize, usize)>> {
    let mut num_add_order = vec![(0, 0); st.ev.n_ev];
    let charge_margin = 1000;
    for eid in 0..st.ev.n_ev {
        perm[eid].clear();
        perm[eid].push((st.ev.next[eid].0, INF));
        perm[eid].extend(st.ev.order[eid].iter().filter(|&&oid| st.order.state[oid] == 1).map(|&oid| (st.order.term_pos[oid].1, oid)));
        if st.ev.charge[eid] < charge_margin + st.ev.delta_ev_move * (st.graph.nearest_grid[st.ev.next[eid].0].1 + st.ev.next[eid].1) {num_add_order[eid] = (perm[eid].len(), 0);}
        else {num_add_order[eid] = (perm[eid].len(), st.ev.n_trans_max + 1 - perm[eid].len());}
    }

    let mut update = true;
    let mut used_order = vec![false; 1024];
    while update {
        update = false;
        let mut mn = INF;//10?
        let mut neid = 0;
        let mut noid = 0;
        for &oid in st.order.ids.iter() {
            if used_order[oid] || st.order.state[oid] == 1 {continue;}
            for eid in 0..st.ev.n_ev {
                if num_add_order[eid].0 + num_add_order[eid].1 * 2 == perm[eid].len() {continue;}
                let pos = match perm[eid].last() {
                    Some(&(v, _o)) => {v},
                    None => {st.ev.next[eid].0}
                };
                if st.graph.costs[st.order.term_pos[oid].0][pos] + st.graph.costs[st.order.term_pos[oid].0][st.order.term_pos[oid].1] < mn {
                    mn = st.graph.costs[st.order.term_pos[oid].0][pos] + st.graph.costs[st.order.term_pos[oid].0][st.order.term_pos[oid].1];
                    neid = eid;
                    noid = oid;
                }
            }
            
        }
        if mn != INF {
            update = true;
            used_order[noid] = true;
            perm[neid].push((st.order.term_pos[noid].0, noid));
            perm[neid].push((st.order.term_pos[noid].1, noid));
        }
    }
    for eid in 0..st.ev.n_ev {
        num_add_order[eid].1 = (perm[eid].len() - num_add_order[eid].0) / 2;
    }
    let mut dq = vec![VecDeque::new(); st.ev.n_ev];
    *perm = optimize_ev_route(st, perm, &num_add_order);
    for eid in 0..st.ev.n_ev {
        for &nx in perm[eid].iter().skip(1) {
            if nx.0 != V && nx.1 != INF {dq[eid].push_back(nx);}
        }
    }
    dq
}

#[inline]
fn move_to(st: &mut State, eid: usize, fr: usize, to: usize, curcharge: &mut usize) {
    let mut pos = fr;
    while pos != to {
        let subnx = st.graph.next[pos][to];
        let dist = st.graph.costs[pos][subnx];
        for _ in 0..dist {
            st.ev.commands[eid].push_back(Command::Move(subnx));
            *curcharge -= st.ev.delta_ev_move;
        }
        pos = subnx;
    }
}
#[inline]
fn charge_from_grid(st: &mut State, eid: usize, gid: usize, target_charge: usize, curcharge: &mut usize) {
    let charge_turn = ((target_charge - *curcharge) + st.ev.v_ev_max - 1) / st.ev.v_ev_max;
    for _ in 0..charge_turn {
        let vc = std::cmp::min(st.ev.v_ev_max, target_charge - *curcharge);
        st.ev.commands[eid].push_back(Command::ChargeFromGrid(vc));
        if st.nano_grid.charge[gid] >= vc {st.nano_grid.charge[gid] -= vc;}
        else {st.nano_grid.charge[gid] = 0;}
        *curcharge += vc;
    }
}
#[inline]
fn charge_to_grid(st: &mut State, eid: usize, gid: usize, target_charge: usize, curcharge: &mut usize){
    let charge_turn = ((*curcharge - target_charge) + st.ev.v_ev_max - 1) / st.ev.v_ev_max;
    for _ in 0..charge_turn {
        let vc = std::cmp::min(st.ev.v_ev_max, *curcharge - target_charge);
        st.ev.commands[eid].push_back(Command::ChargeToGrid(vc));
        st.nano_grid.charge[gid] += vc;
        *curcharge -= vc;
    }
}

fn construct_charge_command_all(st: &mut State, available: &mut Vec<usize>, match_vec: &Vec<usize>, t: usize) {
    let mut cnt = vec![0; V];
    for &mt in match_vec.iter() {
        if mt != V {cnt[mt] += 1;}
    }
    for eid in 0..st.ev.n_ev {
        st.ev.commands[eid].clear();
        if match_vec[eid] == V {continue;}
        let mut curcharge = st.ev.charge[eid];
        let pos = st.ev.next[eid].0; 
        for _ in 0..st.ev.next[eid].1 {
            st.ev.commands[eid].push_back(Command::Move(pos));
            curcharge -= st.ev.delta_ev_move;
        }
        move_to(st, eid, pos, match_vec[eid], &mut curcharge);
        let gid = st.nano_grid.pos_to_id[match_vec[eid]];
        let target_charge = std::cmp::min(25000, curcharge + available[gid] / cnt[match_vec[eid]]);
        //target_charge = target_charge.min(st.ev.c_ev_max - curcharge);
        let amount = std::cmp::min(available[gid], target_charge - curcharge);
        available[gid] -= amount;
        st.nano_grid.overflow[gid] -= (target_charge - curcharge) as i32;
        charge_from_grid(st, eid, gid, target_charge, &mut curcharge);
        st.nano_grid.chrg_end_time[gid].push_back(t + st.ev.commands[eid].len());
    }
}

fn ev_charge_command (st: &mut State, available: &mut Vec<usize>, eid: usize, curcharge: &mut usize, nx_pos: usize, que_len: usize, t: usize, mini: usize) -> usize {
    let pos = st.ev.next[eid].0;
    let ini_cost = st.ev.next[eid].1;
    let mut mn = mini;
    let mgn = 1000;
    let mut ngid = V;
    for gid in 0..st.nano_grid.n_grid {
        let tav = std::cmp::max(available[gid], st.nano_grid.overflow[gid] as usize);
        if tav < mgn {continue;}
        let subnx = st.nano_grid.x[gid];
        let subcost = st.graph.costs[pos][subnx];
        if st.nano_grid.chrg_end_time[gid].len() == 2 {
            let &fr = st.nano_grid.chrg_end_time[gid].front().unwrap();
            let &bk = st.nano_grid.chrg_end_time[gid].back().unwrap();
            if t + subcost + ini_cost < std::cmp::min(fr, bk) {continue;}
        }
        if (subcost + ini_cost) * st.ev.delta_ev_move > *curcharge || tav < (subcost + ini_cost) * st.ev.delta_ev_move + 500 {continue;}
        let diff_cost = subcost + st.graph.costs[subnx][nx_pos] - st.graph.costs[pos][nx_pos];
        
        if diff_cost < mini && st.nano_grid.overflow[gid] > 1000 {
            ngid = gid;
            break;
        }
        
        if diff_cost < mn {
            mn = diff_cost;
            ngid = gid;
        }   
    }
    if ngid == V {return V;}
    let limit = if que_len <= 4 || st.nano_grid.in_sudden[ngid] || st.nano_grid.overflow[ngid] > 1000 {st.ev.c_ev_max} else {3000};
    //let limit = st.ev.c_ev_max;
    for _ in 0..ini_cost {
        st.ev.commands[eid].push_back(Command::Move(pos));
        *curcharge -= st.ev.delta_ev_move;
    }
    move_to(st, eid, pos, st.nano_grid.x[ngid], curcharge);
    let av = limit.min(std::cmp::max(available[ngid], st.nano_grid.overflow[ngid].max(0) as usize));
    //let av = limit.min(available[ngid]);
    let target_charge = std::cmp::min(*curcharge + av, *curcharge + (24600 - *curcharge + st.ev.v_ev_max - 1) / st.ev.v_ev_max * st.ev.v_ev_max);
    let amount = std::cmp::min(available[ngid], target_charge - *curcharge);
    available[ngid] -= amount;
    st.nano_grid.overflow[ngid] -= (target_charge - *curcharge) as i32;
    charge_from_grid(st, eid, ngid, target_charge, curcharge);

    if st.nano_grid.chrg_end_time[ngid].len() == 2 { 
        let &fr = st.nano_grid.chrg_end_time[ngid].front().unwrap();
        let &bk = st.nano_grid.chrg_end_time[ngid].back().unwrap();
        if fr < bk {
            st.nano_grid.chrg_end_time[ngid].pop_front();
        }
        else {
            st.nano_grid.chrg_end_time[ngid].pop_back();
        }
    }
    st.nano_grid.chrg_end_time[ngid].push_back(t + st.ev.commands[eid].len());

    ngid
}

fn pick_order_command (st: &mut State, eid: usize, curcharge: &mut usize, nx_pos: usize, _t: usize, mini: usize) -> usize {
    let mut noid = INF;
    let mut mn = mini;
    let pos = st.ev.next[eid].0;
    let ini_cost = st.ev.next[eid].1;
    for &oid in st.order.ids.iter() {
        if st.order.assign_ev[oid] != V {continue;}
        let nx = st.order.term_pos[oid].0;
        let diff_cost = st.graph.costs[pos][nx] + st.graph.costs[nx][nx_pos] - st.graph.costs[pos][nx_pos];
        if mn > diff_cost && (ini_cost + st.graph.costs[pos][nx]) * st.ev.delta_ev_move + 500 <= *curcharge {
            mn = diff_cost;
            noid = oid;
        }
    }
    if noid != INF {
        st.order.assign_ev[noid] = eid;
        for _ in 0..ini_cost {
            st.ev.commands[eid].push_back(Command::Move(pos));
            *curcharge -= st.ev.delta_ev_move;
        }
        move_to(st, eid, pos, st.order.term_pos[noid].0, curcharge);
        st.ev.commands[eid].push_back(Command::Pickup(noid));
    }
    
    noid
    
}

fn grid_charge_command (st: &mut State, gid: usize, eid: usize, curcharge: &mut usize, t: usize) -> usize {
    let pos = st.ev.next[eid].0;
    let ini_cost = st.ev.next[eid].1;
    let ngid = gid;
    let subnx = st.nano_grid.x[gid];
    let subcost = st.graph.costs[pos][subnx];
    if st.nano_grid.chrg_end_time[gid].len() == 2 {
        let &fr = st.nano_grid.chrg_end_time[gid].front().unwrap();
        let &bk = st.nano_grid.chrg_end_time[gid].back().unwrap();
        if t + subcost + ini_cost < std::cmp::min(fr, bk) {return V;}
    }
    let mgn = if st.nano_grid.in_sudden[gid] {15000} else {5000};
    if (subcost + ini_cost) * st.ev.delta_ev_move + mgn > *curcharge {return V;}

    for _ in 0..ini_cost {
        st.ev.commands[eid].push_back(Command::Move(pos));
        *curcharge -= st.ev.delta_ev_move;
    }
    move_to(st, eid, pos, st.nano_grid.x[ngid], curcharge);

    let nd = std::cmp::min(*curcharge - 1500, st.nano_grid.need[ngid] as usize);
    let target_charge = std::cmp::max(*curcharge - nd, (*curcharge - 1500 + st.ev.v_ev_max - 1) / st.ev.v_ev_max * st.ev.v_ev_max);
    st.nano_grid.need[ngid] -= (*curcharge - target_charge) as i32;
    charge_to_grid(st, eid, ngid, target_charge, curcharge);

    if st.nano_grid.chrg_end_time[gid].len() == 2 { 
        let &fr = st.nano_grid.chrg_end_time[gid].front().unwrap();
        let &bk = st.nano_grid.chrg_end_time[gid].back().unwrap();
        if fr < bk {
            st.nano_grid.chrg_end_time[gid].pop_front();
        }
        else {
            st.nano_grid.chrg_end_time[gid].pop_back();
        }
    }
    st.nano_grid.chrg_end_time[gid].push_back(t + st.ev.commands[eid].len());

    ngid
}

fn avoid_overflow_command (st: &mut State, available: &mut Vec<usize>, gid: usize, eid: usize, curcharge: &mut usize, t: usize) -> usize {
    let pos = st.ev.next[eid].0;
    let ini_cost = st.ev.next[eid].1;
    let ngid = gid;
    let subnx = st.nano_grid.x[gid];
    let subcost = st.graph.costs[pos][subnx];
    if st.nano_grid.chrg_end_time[gid].len() == 2 {
        let &fr = st.nano_grid.chrg_end_time[gid].front().unwrap();
        let &bk = st.nano_grid.chrg_end_time[gid].back().unwrap();
        if t + subcost + ini_cost < std::cmp::min(fr, bk) {return V;}
    }
    if (subcost + ini_cost) * st.ev.delta_ev_move > *curcharge {return V;}

    for _ in 0..ini_cost {
        st.ev.commands[eid].push_back(Command::Move(pos));
        *curcharge -= st.ev.delta_ev_move;
    }
    move_to(st, eid, pos, st.nano_grid.x[ngid], curcharge);
    let av = std::cmp::max(available[ngid], st.nano_grid.overflow[ngid].max(0) as usize);
    //let av = std::cmp::min(available[ngid], st.nano_grid.overflow[ngid].max(0) as usize);
    let target_charge = std::cmp::min(*curcharge + av, *curcharge + (24600 - *curcharge + st.ev.v_ev_max - 1) / st.ev.v_ev_max * st.ev.v_ev_max);
    let amount = std::cmp::min(available[ngid], target_charge - *curcharge);
    available[ngid] -= amount;
    st.nano_grid.overflow[ngid] -= (target_charge - *curcharge) as i32;
    charge_from_grid(st, eid, ngid, target_charge, curcharge);

    if st.nano_grid.chrg_end_time[gid].len() == 2 { 
        let &fr = st.nano_grid.chrg_end_time[gid].front().unwrap();
        let &bk = st.nano_grid.chrg_end_time[gid].back().unwrap();
        if fr < bk {
            st.nano_grid.chrg_end_time[gid].pop_front();
        }
        else {
            st.nano_grid.chrg_end_time[gid].pop_back();
        }
    }
    st.nano_grid.chrg_end_time[gid].push_back(t + st.ev.commands[eid].len());

    ngid
}

fn main () {
    let timer = Timer::new();

    let n_solution = read!(usize);
    //graph inputs
    let (v, e) = read!(usize, usize);
    let mut gr = Graph {v: v, e: e, graph: vec![Vec::new(); v], costs: vec![vec![INF; v+1]; v+1], next:vec![vec![INF; v]; v], nearest_grid: vec![(0, 0); v]};
    
    for _ in 0..e {
        let (ui, vi, di) = read!(usize1, usize1, usize);
        gr.graph[ui].push(vi);
        gr.graph[vi].push(ui);
        gr.costs[ui][vi] = di;
        gr.costs[vi][ui] = di;
    }

    for i in 0..v {
        gr.costs[i][i] = 0;
        gr.costs[i][v] = 0;
        gr.costs[v][i] = 0;
        for j in 0..v {
            gr.next[i][j] = j;
        }
    }
    for k in 0..v {
        for i in 0..v {
            for j in 0..v {
                if gr.costs[i][j] > gr.costs[i][k] + gr.costs[k][j] {
                    gr.costs[i][j] = gr.costs[i][k] + gr.costs[k][j];
                    gr.next[i][j] = gr.next[i][k];
                }
            }
        }
    }

    //nanogrid inputs
    let daytype = match read!(usize) {
        0 => DayType::Fine,
        1 => DayType::FineSquall,
        2 => DayType::Rainy,
        _ => DayType::RainyFine,
    };

    let (n_div, n_pattern, sigma2, p_event, d_event) = read!(usize, usize, usize, f64, i32);
    let pw_pattern = read!([i32]; n_pattern);
    let (n_grid, c_init, c_max, v_max) = read!(usize, usize, usize, usize);
    let mut pw_expected = vec![vec![0; n_div+1]; n_pattern];
    for i in 0..n_pattern {
        for j in 1..n_div+1 {
            pw_expected[i][j] = pw_expected[i][j-1] + pw_pattern[i][j-1] * DIV_RANGE as i32;
        }
    }
    let mut fsthalf_minexpect_divt = vec![0; n_pattern];
    let mut maxexpect_divt = vec![0; n_pattern];
    for i in 0..n_pattern {
        let mut mnid = 0;
        let mut mxid = 0;
        for j in 1..n_div+1 {
            if j <= 10 && pw_expected[i][mnid] > pw_expected[i][j] {mnid = j;}
            //if pw_expected[i][mnid] > pw_expected[i][j] {mnid = j;}
            if pw_expected[i][mxid] < pw_expected[i][j] {mxid = j;}
        }
        fsthalf_minexpect_divt[i] = mnid;
        maxexpect_divt[i] = mxid;
    }
    let mut nano_grid = NanoGrids {
        n_grid: n_grid, n_div: n_div, p_event: p_event, delta_event: d_event, n_pattern: n_pattern,
        pw_predicts: pw_pattern, sigma_ele: sigma2,
        c_grid_init: c_init, c_grid_max: c_max, v_grid_max: v_max,
        x: Vec::new(),
        patterns: Vec::new(),
        day_type: daytype,
        charge: vec![c_init; n_grid],
        pos_to_id: vec![v; v],
        pw_expected: pw_expected,
        visit_ev_id: vec![VecDeque::new(); n_grid],
        chrg_end_time: vec![VecDeque::new(); n_grid],
        use_count: vec![0; n_grid],
        in_sudden: vec![false; n_grid],
        overflow: vec![0; n_grid],
        need: vec![0; n_grid],
        fsthalf_minexpect_divt: fsthalf_minexpect_divt,
        maxexpect_divt: maxexpect_divt,
    };
    for i in 0..n_grid {
        let (xi, pi) = read!(usize, usize);
        nano_grid.x.push(xi - 1);
        nano_grid.patterns.push(pi - 1);
        nano_grid.pos_to_id[xi - 1] = i;
    }
    //nearest grid
    for i in 0..V {
        let mut mn = INF;
        let mut nid = 0;
        for &pos in nano_grid.x.iter() {
            if gr.costs[i][pos] < mn {
                mn = gr.costs[i][pos];
                nid = pos;
            }
        }
        gr.nearest_grid[i] = (nid, mn);
    }

    //EV inputs
    let (n_ev, c_init, c_max, v_max, n_trans_max, delta_move) = read!(usize, usize, usize, usize, usize, usize);
    let ev_pos = read!(usize1; n_ev);
    let next: Vec<(usize, usize)> = ev_pos.iter().map(|&x| (x, 0usize)).collect();
    let ev = EVs {
        n_ev: n_ev,
        c_ev_init: c_init, c_ev_max: c_max, v_ev_max: v_max, n_trans_max: n_trans_max, delta_ev_move: delta_move,
        charge: vec![c_init; n_ev],
        pos: ev_pos,
        next: next,
        adj: vec![Vec::new(); n_ev],
        commands: vec![VecDeque::new(); n_ev],
        pre_commands: vec![Command::Stay; n_ev],
        order: vec![Vec::new(); n_ev],
        order_cnt: vec![0; n_ev],
        to_charge_flag: vec![false; n_ev],
    };

    let order = Orders {
        n_order : 0,
        ids: Vec::new(),
        term_pos: vec![(0, 0); 1024],
        state: vec![0; 1024],
        time: vec![0; 1024],
        assign_ev: vec![V; 1024],
    };

    //other inputs
    let (_p_const_trans, _t_last) = read!(f64, usize);
    let (_penalty, _gamma, _score_ref_ele, _score_ref_trans) = read!(usize, f64, f64, f64);
    let _tmax = read!(usize);

    //init
    let mut ln16: [f64; 65536] = [0.0; 65536];
    for i in 0..65536 {
        ln16[i] = (i as f64 / 65536.0 + 1.0/(2.0*65536.0)).ln();
    }
    let mut state = State
    {
        rnd: XorShift32::new(0),
        timer: timer,
        ln: ln16,
        n_solution: n_solution,
        graph: gr,
        nano_grid: nano_grid,
        ev: ev,
        order: order,
        cur_solution_id: 0,
        sudden: false,
    };

    let mut scorelist = Vec::new();
    let stdout = std::io::stdout();
    let mut stdout = BufWriter::new(stdout.lock());
    let mut available = vec![0; state.nano_grid.n_grid];
    let mut perm = vec![Vec::new(); state.ev.n_ev];
    
    let check_t = match state.nano_grid.day_type {
        DayType::Fine | DayType::Rainy => {0},
        DayType::FineSquall | DayType::RainyFine => {1}
    };
    
    let terms = match state.nano_grid.day_type {
        DayType::Fine | DayType::FineSquall => {
            match state.ev.n_ev {
                0..=16 => {vec![100,100,40,100,40]},
                17..=19 => {vec![100,100,50,100,50]},
                20..=22 => {vec![100,100,100,100,60]},
                _ => {vec![100,100,100,100,75]}
            }
        },
        DayType::Rainy | DayType::RainyFine => {
            match state.ev.n_ev {
                0..=16 => {vec![100,100,100,100,40]},
                17..=22 => {vec![100,100,100,100,60]},
                _ => {vec![100,100,100,100,75]}
            }
        }
    };
    
    //let mut additional_pick = vec![0usize; state.ev.n_ev];
    //let term = 100;
    for n in 0..state.n_solution {
        state.cur_solution_id = n;
        if n > 0 {init_state(&mut state);}
        let term = terms[n];
        state.sudden = false;
        let mut stopflag = false;

        let mut pdq = vec![VecDeque::new(); state.ev.n_ev];
        let mut nx_calc = check_t;
        for t in 0..T_MAX {
            read_info(&mut state, t / DIV_RANGE);
            //if OUT_INFO && t % DIV_RANGE == 0 {print_status(&state, t);}
            
            //if t % DIV_RANGE == check_t {
            if t % DIV_RANGE == check_t {
                calc_need_charge(&mut state, t);
            }

            //if t % term == check_t {
            if t == nx_calc {//|| state.sudden && t % DIV_RANGE == check_t {
                nx_calc = t + term;
                stopflag = false;
                calc_expected_charge(&mut state, &mut available, t);
                let bd = 1;
                let bl = match state.nano_grid.day_type {
                    DayType::Fine | DayType::FineSquall => {17},
                    DayType::Rainy | DayType::RainyFine => {15}
                };
                if t == check_t {
                    let match_vec = optimaize_init_matching(&mut state, &available);
                    construct_charge_command_all(&mut state, &mut available, &match_vec, t);
                    for eid in 0..state.ev.n_ev {
                        state.ev.to_charge_flag[eid] = true;
                    }
                }
                else if n <= bd {
                    let mut flag = true;
                    match n {
                        0 => {
                            for pat in 0..state.nano_grid.n_pattern {
                                if state.nano_grid.fsthalf_minexpect_divt[pat] > t / DIV_RANGE {flag = false;}
                            }
                            let mut av_sum = 0;
                            for i in 0..state.nano_grid.n_grid {
                                //if available[i] > 1000 {flag = false;}
                                av_sum += available[i];
                            }
                            if av_sum > 5000 {flag = false;}
                        }
                        1 => {
                            for _pat in 0..state.nano_grid.n_pattern {
                                //if state.nano_grid.maxexpect_divt[pat] > t / DIV_RANGE {flag = false;}
                                if bl > t / DIV_RANGE {flag = false;}
                            }
                        },
                        _ => {}
                    }
                    
                    if flag {
                        stopflag = true;
                    }
                    else {
                        for ev in 0..state.ev.n_ev {
                            state.ev.commands[ev].clear();
                        }
                        pdq = construct_ev_route(&mut state, &mut perm, t);
                    }
                }
                else {
                    for ev in 0..state.ev.n_ev {
                        state.ev.commands[ev].clear();
                    }
                    pdq = construct_ev_route(&mut state, &mut perm, t);
                }
                if !stopflag {
                    
                    for gid in 0..state.nano_grid.n_grid {
                        state.nano_grid.chrg_end_time[gid].clear();
                    }
                    
                    for i in state.order.assign_ev.iter_mut() {*i = V;}
                    for i in state.ev.order_cnt.iter_mut() {*i = 0;}
                    //for i in additional_pick.iter_mut() {*i = 0;}
                    for (eid, ev) in perm.iter().enumerate() {
                        for &(pos, oid) in ev.iter() {
                            if oid > 1024 {continue;}
                            if state.order.term_pos[oid].1 == pos {state.ev.order_cnt[eid] += 1;}
                            state.order.assign_ev[oid] = eid; 
                        }
                    }
                }
            }

            for eid in 0..state.ev.n_ev {
                if state.ev.commands[eid].is_empty() && t > check_t {
                    state.ev.to_charge_flag[eid] = true;
                    if let Some(&nx) = pdq[eid].front() {
                        let curpos = state.ev.next[eid].0;
                        let mut curcharge = state.ev.charge[eid];
                        let mut smpos_flag = false;
                        if curpos == nx.0 {
                            smpos_flag = true;
                        }
                        
                        if let DayType::FineSquall = state.nano_grid.day_type {
                            let mut tnx = V;
                            if !smpos_flag && curcharge > 15000 && pdq[eid].len() <= 3 {
                                for gid in 0..state.nano_grid.n_grid {
                                    let gpos = state.nano_grid.x[gid];
                                    let diff_cost = state.graph.costs[curpos][gpos] + state.graph.costs[gpos][nx.0] - state.graph.costs[curpos][nx.0];
                                    if diff_cost <= 5 && state.nano_grid.in_sudden[gid] && state.nano_grid.need[gid] > 1000 {
                                        tnx = grid_charge_command(&mut state, gid, eid, &mut curcharge, t);
                                        if tnx != V {break;}
                                    }
                                }
                                if tnx != V {
                                    continue;
                                }
                            }
                        }
                        
                        if let DayType::RainyFine = state.nano_grid.day_type {
                            let mut tnx = V;
                            if !smpos_flag && curcharge < 20000 && pdq[eid].len() <= 3 {
                                for gid in 0..state.nano_grid.n_grid {
                                    let gpos = state.nano_grid.x[gid];
                                    let diff_cost = state.graph.costs[curpos][gpos] + state.graph.costs[gpos][nx.0] - state.graph.costs[curpos][nx.0];
                                    if diff_cost <= 10 && state.nano_grid.in_sudden[gid] && state.nano_grid.overflow[gid] > 3000 {
                                        tnx = avoid_overflow_command(&mut state, &mut available, gid, eid, &mut curcharge, t);
                                        if tnx != V {break;}
                                    }
                                }
                                if tnx != V {
                                    continue;
                                }
                            }
                        }

                        
                        if !smpos_flag && curcharge > 5000 && pdq[eid].len() <= 3 {
                            let mut tnx = V;
                            for gid in 0..state.nano_grid.n_grid {
                                let gpos = state.nano_grid.x[gid];
                                let diff_cost = state.graph.costs[curpos][gpos] + state.graph.costs[gpos][nx.0] - state.graph.costs[curpos][nx.0];
                                if diff_cost <= 5 && state.nano_grid.need[gid] > 1000 {
                                    tnx = grid_charge_command(&mut state, gid, eid, &mut curcharge, t);
                                    if tnx != V {break;}
                                }
                            }
                            if tnx != V {
                                continue;
                            }
                        }

                        if !smpos_flag && t < T_MAX - 100 && state.ev.charge[eid] < 1500 {//quelen次第でcharge制限かける
                            let tnx = ev_charge_command(&mut state, &mut available, eid, &mut curcharge, nx.0, pdq[eid].len(), t, 15);
                            if tnx != V {
                                continue;
                            }
                        }
                        /*
                        if !stopflag && state.ev.order_cnt[eid] < state.ev.n_trans_max && additional_pick[eid] == 0 {
                            let noid = pick_order_command(&mut state, eid, &mut curcharge, nx.0, t, 1);
                            if noid != INF {
                                state.ev.order_cnt[eid] += 1;
                                additional_pick[eid] += 1;
                                //insert dst
                                let mut pre = state.order.term_pos[noid].0;
                                let dst = state.order.term_pos[noid].1;
                                let mut ins_id = INF;
                                let mut mn = INF;
                                for (i, &(v, _)) in pdq[eid].iter().enumerate() {
                                    let diff = state.graph.costs[pre][dst] + state.graph.costs[dst][v] - state.graph.costs[pre][v];
                                    if diff < mn {
                                        mn = diff;
                                        ins_id = i;
                                    }
                                    pre = v; 
                                }
                                if mn < 3 {
                                    pdq[eid].insert(ins_id, (dst, noid));
                                }
                                else {
                                    pdq[eid].push_back((dst, noid));
                                }
                                continue;
                            }
                            
                        }
                        */
                        //charge >24600 
                        if !smpos_flag && t < T_MAX - 100 && state.ev.charge[eid] < 24600 {// && pdq[eid].len() <= 3
                            let tnx = ev_charge_command(&mut state, &mut available, eid, &mut curcharge, nx.0, pdq[eid].len(), t, 5);
                            if tnx != V {
                                continue;
                            }
                        }

                        //なければ次のnx.0
                        state.ev.to_charge_flag[eid] = false;

                        let pos = state.ev.next[eid].0;
                        let ini_cost = state.ev.next[eid].1;
                        if (ini_cost + state.graph.costs[pos][nx.0]) * state.ev.delta_ev_move > curcharge {continue;}
                        for _ in 0..ini_cost {
                            state.ev.commands[eid].push_back(Command::Move(state.ev.next[eid].0));
                            curcharge -= state.ev.delta_ev_move;
                        }
                        move_to(&mut state, eid, pos, nx.0, & mut curcharge);
                        if nx.1 < 1024 && state.order.term_pos[nx.1].0 == nx.0 {
                            state.ev.commands[eid].push_back(Command::Pickup(nx.1));
                            pdq[eid].pop_front();
                        }
                        else {
                            let mut cnt = 0;
                            for &tnx in pdq[eid].iter() {
                                if nx.0 == tnx.0 && tnx.1 < 1024 && state.order.term_pos[tnx.1].1 == nx.0 {
                                    cnt += 1;
                                }
                                else {break;}
                            }
                            for _ in 0..cnt {
                                state.ev.order_cnt[eid] -= 1;
                                pdq[eid].pop_front();
                            }
                        }

                    }
                    else {
                        //que empty
                        let curpos = state.ev.next[eid].0;
                        let mut curcharge = state.ev.charge[eid];

                        if let DayType::FineSquall = state.nano_grid.day_type {
                            let mut tnx = V;
                            if curcharge > 15000 {
                                for gid in 0..state.nano_grid.n_grid {
                                    let gpos = state.nano_grid.x[gid];
                                    let diff_cost = state.graph.costs[curpos][gpos];
                                    if diff_cost <= 15 && state.nano_grid.in_sudden[gid] && state.nano_grid.need[gid] > 1000 {
                                        tnx = grid_charge_command(&mut state, gid, eid, &mut curcharge, t);
                                        if tnx != V {break;}
                                    }
                                }
                                if tnx != V {
                                    continue;
                                }
                            }
                        }
                        
                        if let DayType::RainyFine = state.nano_grid.day_type {
                            let mut tnx = V;
                            if curcharge < 22000 {
                                for gid in 0..state.nano_grid.n_grid {
                                    let gpos = state.nano_grid.x[gid];
                                    let diff_cost = state.graph.costs[curpos][gpos];
                                    if diff_cost <= 10 && state.nano_grid.in_sudden[gid] && state.nano_grid.overflow[gid] > 3000 {
                                        tnx = avoid_overflow_command(&mut state, &mut available, gid, eid, &mut curcharge, t);
                                        if tnx != V {break;}
                                    }
                                }
                                if tnx != V {
                                    continue;
                                }
                            }
                        }

                        if curcharge > 3000 {
                            let mut tnx = V;
                            for gid in 0..state.nano_grid.n_grid {
                                let gpos = state.nano_grid.x[gid];
                                let diff_cost = state.graph.costs[curpos][gpos];
                                if diff_cost <= 15 && state.nano_grid.need[gid] > 1000 {
                                    tnx = grid_charge_command(&mut state, gid, eid, &mut curcharge, t);
                                    if tnx != V {break;}
                                }
                            }
                            if tnx != V {
                                continue;
                            }
                        }

                        if t < T_MAX - 100 && state.ev.charge[eid] < 1500 {
                            let tnx = ev_charge_command(&mut state, &mut available, eid, &mut curcharge, V, 0, t, 15);
                            if tnx != V {
                                continue;
                            }
                        }
                        /*
                        if !stopflag && state.ev.order_cnt[eid] < state.ev.n_trans_max && additional_pick[eid] == 0 {
                            let noid = pick_order_command(&mut state, eid, &mut curcharge, V, t, 3);
                            if noid != INF {
                                state.ev.order_cnt[eid] += 1;
                                additional_pick[eid] += 1;
                                //insert dst
                                let mut pre = state.order.term_pos[noid].0;
                                let dst = state.order.term_pos[noid].1;
                                let mut ins_id = INF;
                                let mut mn = INF;
                                for (i, &(v, _)) in pdq[eid].iter().enumerate() {
                                    let diff = state.graph.costs[pre][dst] + state.graph.costs[dst][v] - state.graph.costs[pre][v];
                                    if diff < mn {
                                        mn = diff;
                                        ins_id = i;
                                    }
                                    pre = v; 
                                }
                                if mn < 5 {
                                    pdq[eid].insert(ins_id, (dst, noid));
                                }
                                else {
                                    pdq[eid].push_back((dst, noid));
                                }
                                continue;
                            }
                        }
                        */
                        //charge >24600 
                        if t < T_MAX - 100 && state.ev.charge[eid] < 24600 {
                            let tnx = ev_charge_command(&mut state, &mut available, eid, &mut curcharge, V, 0, t, 5);
                            if tnx != V {
                                continue;
                            }
                        }

                        state.ev.to_charge_flag[eid] = false;
                    }
                }
            }
            
            let tt = t;
            for dq in state.nano_grid.chrg_end_time.iter_mut() {
                if let Some(&ft) = dq.front() {if tt+1 == ft {dq.pop_front();}}
                if let Some(&ft) = dq.back() {if tt+1 == ft {dq.pop_back();}}
            }
            
            for (ev, com) in state.ev.commands.iter_mut().enumerate() {
                match com.front() {
                    Some(c) => {
                        output_command(c, &mut stdout);
                        state.ev.pre_commands[ev] = c.clone();
                        com.pop_front();
                    },
                    None => {
                        output_command(&Command::Stay, &mut stdout);
                        state.ev.pre_commands[ev] = Command::Stay.clone();
                    }
                }
                
                //print!("stay\n");
            }
            stdout.flush().unwrap();
        }
        //finally :Tmax info
        read_info(&mut state, 19);
        if OUT_INFO {print_status(&state, 1000);}
        let judge_score = read!(f64, f64);
        eprint!("num_rem_order: {}\n", state.order.n_order);
        //if OUT_INFO {eprintln!("{} {}", judge_score.0 , judge_score.1);}
        scorelist.push(judge_score);
    }

    

    
    eprint!("{}\n", state.timer.get_time());
    for i in 0..n_solution {
        eprint!("{} ", scorelist[i].0);
    }
    eprint!("\n");
    for i in 0..n_solution {
        eprint!("{} ", scorelist[i].1);
    }
    /*
    let total_score = read!(f64);
    eprint!("\n{}\n", total_score);
    */
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

