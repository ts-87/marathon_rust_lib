#![allow(unused_imports)]
#![allow(dead_code)]

use std::{collections::VecDeque, io::prelude::*, mem::swap};

const OUT_INFO: bool = true;
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
const INF: usize = std::usize::MAX;
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
    fsthalf_minexpect_divt: Vec<usize>,
    maxexpect_divt: Vec<usize>
}

struct EVs {
    n_ev: usize,
    c_ev_init: usize, c_ev_max: usize, v_ev_max: usize,
    delta_ev_move: usize,
    charge: Vec<usize>,
    pos: Vec<usize>,
    next: Vec<(usize, usize)>,//id,dist
    adj: Vec<Vec<usize>>,
    commands: Vec<VecDeque<Command>>,
    pre_commands: Vec<Command>,
}
#[derive(Clone)]
enum Command {
    Stay,
    Move(usize),
    ChargeFromGrid(usize),
    ChargeToGrid(usize)
}

struct State {
    rnd: XorShift32,
    timer: Timer,
    ln: [f64; 65536],
    graph: Graph,
    nano_grid: NanoGrids,
    ev: EVs,
    sudden: bool,
    complete: bool,
    unblanced: bool,
    after_sudden: bool,
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
    eprint!("\n");
}

fn read_info(st: &mut State, div_t: usize) {
    st.sudden = false;
    let is_sudden_daytype = match st.nano_grid.day_type {
        DayType::FineSquall | DayType::RainyFine => {true},
        _ => {false}
    };
    for i in 0..st.nano_grid.n_grid {
        let (xi, yi, pw_actual, _pw_exceed, _pw_buy) = read!(usize1, usize, i32, usize, usize);
        let pat = st.nano_grid.patterns[i];
        if is_sudden_daytype && (st.nano_grid.pw_predicts[pat][div_t] - pw_actual).abs() > 600 {
            st.sudden = true;
            st.nano_grid.in_sudden[i] = true;
            st.after_sudden = true;
        }
        else {
            st.nano_grid.in_sudden[i] = false;
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
    }
}

/*
fine
-247,-209,-163,-57,5,119,146,200,208,204,190,131,69,22,-95,-139,-173,-165,-97,-68
8,43,88,129,166,93,82,49,56,53,81,60,32,-2,-42,-166,-209,-197,-174,-150
-156,-116,-115,-89,-52,-26,-1,24,48,60,50,20,-43,-78,-115,-114,-157,-140,-123,-118
sum expected
25000  12650   2200  -5950  -8800  -8550  -2600   4700  14700  25100  35300  44800  51350  54800  55900  51150  44200  35550  27300  22450  19050 
25000  25400  27550  31950  38400  46700  51350  55450  57900  60700  63350  67400  70400  72000  71900  69800  61500  51050  41200  32500  25000 
25000  17200  11400   5650   1200  -1400  -2700  -2750  -1550    850   3850   6350   7350   5200   1300  -4450 -10150 -18000 -25000 -31150 -37050
rainy
27,36,47,40,73,82,51,52,56,54,10,7,-32,-52,-137,-146,-145,-89,-61,-32
-38,-31,-22,-14,-7,16,-13,-58,-57,-57,-60,-64,-29,-36,-44,-53,-61,-41,-18,6
-122,-110,-105,-112,-64,-63,-56,-44,-64,-50,30,16,-25,-32,-37,-41,-50,-24,-47,-42
sum expected
25000  26350  28150  30500  32500  36150  40250  42800  45400  48200  50900  51400  51750  50150  47550  40700  33400  26150  21700  18650  17050 
25000  23100  21550  20450  19750  19400  20200  19550  16650  13800  10950   7950   4750   3300   1500   -700  -3350  -6400  -8450  -9350  -9050 
25000  18900  13400   8150   2550   -650  -3800  -6600  -8800 -12000 -14500 -13000 -12200 -13450 -15050 -16900 -18950 -21450 -22650 -25000 -27100 
*/

fn calc_need_priority_fine (st:&mut State, div_t: usize, check_t: usize) -> Vec<(i32, usize)> {
    let n_grid = st.nano_grid.n_grid;
    let mut need = vec![(0i32, 0usize); n_grid];
    st.complete = false;
    let mut tot_ev_charge = 0;
    let mut tot_need_sum = 0;
    for &i in st.ev.charge.iter() {
        tot_ev_charge += i;
    }
    let mut posneedsum = 0;
    let mut fcnt = 0;
    for i in 0..n_grid {
        let mut tmp = -(st.nano_grid.charge[i] as i32);
        let pat = st.nano_grid.patterns[i];
        if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
            tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]];
            //if st.nano_grid.fsthalf_minexpect_divt[pat] == st.nano_grid.n_div {fcnt += 1;}
        }
        else {
            tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.n_div];
            fcnt += 1;
        }
        if check_t == 1 {
            tmp += st.nano_grid.pw_predicts[pat][div_t];
            if st.nano_grid.in_sudden[i] {
                tmp += st.nano_grid.delta_event * (DIV_RANGE as i32);
            }
        }

        let add_c = match st.nano_grid.day_type {
            DayType::Fine | DayType::Rainy => {250},
            DayType::FineSquall => {0},
            DayType::RainyFine => {500}
        };
        tmp += add_c;

        if !st.nano_grid.in_sudden[i] {
            let stmn = (st.nano_grid.c_grid_max - st.nano_grid.charge[i]) as i32;
            tmp = std::cmp::min(tmp, stmn);
        }
        need[i] = (tmp, st.nano_grid.x[i]);
        if need[i].0 < 0 {
            need[i].0 = need[i].0.max(-(st.nano_grid.charge[i] as i32) + add_c);
            
        }
        if tmp > 0 {posneedsum += tmp;}
        tot_need_sum += tmp;
    }
    if fcnt == n_grid && posneedsum == 0 {st.complete = true;}
    st.unblanced = false;
    if (tot_ev_charge as i32) < tot_need_sum + 3000 {
        st.unblanced = true;
    }

    need.sort();
    need.reverse();
    need
}

fn calc_need_priority_rainy (st:&mut State, div_t: usize, check_t: usize) -> Vec<(i32, usize)> {
    let n_grid = st.nano_grid.n_grid;
    let mut need = vec![(0i32, 0usize); n_grid];
    st.complete = false;
    let mut tot_ev_charge = 0;
    let mut tot_need_sum = 0;
    for &i in st.ev.charge.iter() {
        tot_ev_charge += i;
    }
    let mut posneedsum = 0;
    let mut fcnt = 0;
    for i in 0..n_grid {
        let mut tmp = -(st.nano_grid.charge[i] as i32);
        let pat = st.nano_grid.patterns[i];
        if check_t == 2 {//2: deadcode
            if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
                tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]];
                
                if st.nano_grid.fsthalf_minexpect_divt[pat] == st.nano_grid.n_div {
                    fcnt += 1;
                }
            }
            else {
                if st.nano_grid.pw_expected[pat][st.nano_grid.n_div] >= 0 {
                    tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.n_div];
                }
                else {
                    
                    if st.nano_grid.pw_expected[pat][st.nano_grid.n_div.min(div_t + 7)] > 0 {
                        tmp += st.nano_grid.pw_expected[pat][div_t];
                    }
                    else {
                        tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.n_div.min(div_t + 7)];
                    }
                    
                }
                if div_t + 7 >= st.nano_grid.n_div {fcnt += 1;}
            }
            tmp += st.nano_grid.pw_predicts[pat][div_t];
            if st.nano_grid.in_sudden[i] {
                tmp -= st.nano_grid.delta_event * (DIV_RANGE as i32);
            }
        }
        else {
            if div_t < st.nano_grid.fsthalf_minexpect_divt[pat] {
                tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.fsthalf_minexpect_divt[pat]];
                if st.nano_grid.fsthalf_minexpect_divt[pat] == st.nano_grid.n_div {
                    fcnt += 1;
                }
            }
            else {
                tmp += st.nano_grid.pw_expected[pat][div_t] - st.nano_grid.pw_expected[pat][st.nano_grid.n_div];
                fcnt += 1;
            }
            if check_t == 1 {
                tmp += st.nano_grid.pw_predicts[pat][div_t];
                if st.nano_grid.in_sudden[i] {
                    //tmp -= st.nano_grid.delta_event * (DIV_RANGE as i32);
                    if tmp < 0 {tmp -= 780 * (DIV_RANGE as i32);}//v_ev_max * 2
                    else {tmp -= st.nano_grid.delta_event * (DIV_RANGE as i32);}
                }
            }
        }
        let add_c = match st.nano_grid.day_type {
            DayType::Fine | DayType::Rainy => {250},
            DayType::FineSquall => {0},
            DayType::RainyFine => {500}
        };
        tmp += add_c;
        let stmn = (st.nano_grid.c_grid_max - st.nano_grid.charge[i]) as i32;
        tmp = std::cmp::min(tmp, stmn);
        need[i] = (tmp, st.nano_grid.x[i]);
        if need[i].0 < 0 && !st.nano_grid.in_sudden[i] {
            need[i].0 = need[i].0.max(-(st.nano_grid.charge[i] as i32) + add_c);
        }
        
        if tmp > 0 {posneedsum += tmp;}
        tot_need_sum += tmp;
        
    }
    //if fcnt == n_grid && posneedsum == 0 {st.complete = true;}

    st.unblanced = false;
    if (tot_ev_charge as i32) < tot_need_sum + 3000 && div_t >= 10 {
        st.unblanced = true;
    }

    need.sort();
    need.reverse();
    need
}

fn optimize_ev_route (st: &mut State, mut perm: Vec<Vec<usize>>, need: &Vec<i32>, div_t: usize) -> Vec<Vec<usize>> {
    let mut bestscore = 0.;
    let mut grid_penalty = vec![0.; st.nano_grid.n_grid];
    let mut used_charge_sum = vec![0; st.nano_grid.n_grid];

    let (need_weight, avail_weight) = match st.nano_grid.day_type {
        DayType::Fine | DayType::FineSquall => {
            if div_t == 0 {(5., 15.)} else {(4., 12.)}
        }
        DayType::Rainy | DayType::RainyFine => {(4., 12.)},
    };
    let rem = if div_t <= 12 {1000} else {0};

    let mut cntgrid = vec![0; st.nano_grid.n_grid];
    for ev in perm.iter() {
        for j in 1..ev.len() {
            if ev[j] == V {continue;}
            cntgrid[st.nano_grid.pos_to_id[ev[j]]] += 1;
        }
    }
    let mut p_pena = vec![vec![0; 4]; st.nano_grid.n_grid];
    let mut pena_list = vec![vec![0; 4]; st.nano_grid.n_grid];
    let lim_t = 3.min(20 - div_t);
    let lim_step = DIV_RANGE+10;
    for gid in 0..st.nano_grid.n_grid {
        let pat = st.nano_grid.patterns[gid];
        let mut psum = 0;
        for j in 0..lim_t {
            let mut tmp = if j == 0 {st.nano_grid.charge[gid] as i32} else {0}; 
            tmp += st.nano_grid.pw_expected[pat][div_t+j+1] - st.nano_grid.pw_expected[pat][div_t+j];
            if j== 0 && st.sudden {
                match st.nano_grid.day_type {
                    DayType::FineSquall => {tmp -= 50000;}
                    DayType::RainyFine => {tmp += 50000;}
                    _ => {} 
                }
            }
            psum += tmp;
            if psum > st.nano_grid.c_grid_max as i32 {
                p_pena[gid][j+1] = p_pena[gid][j] + psum - st.nano_grid.c_grid_max as i32;
                psum = st.nano_grid.c_grid_max as i32;
            }
            else if psum > 0 {
                p_pena[gid][j+1] = p_pena[gid][j];
            }
            else if psum < 0 {
                p_pena[gid][j+1] = p_pena[gid][j] + psum;
                psum = 0;
            }
        }
    }
    for (eid, ev) in perm.iter().enumerate() {
        let mut curcharge = st.ev.charge[eid];
        let mut cost = st.ev.next[eid].1 * st.ev.delta_ev_move;
        curcharge -= cost;
        let mut pre = st.ev.next[eid].0;
        let mut cur_step = st.ev.next[eid].1;
        for j in 1..ev.len() {
            let nxcost = st.graph.costs[pre][ev[j]] * st.ev.delta_ev_move;
            if curcharge < nxcost {break;}
            //cost += nxcost;
            //curcharge -= nxcost;
            if ev[j] != V {
                cost += nxcost;
                curcharge -= nxcost;
                cur_step += st.graph.costs[pre][ev[j]];
                let mut tmp = 0;
                let gid = st.nano_grid.pos_to_id[ev[j]];
                if need[gid] >= 0 {
                    if curcharge >= 1000 {
                        tmp = std::cmp::min(curcharge as i32 - rem, need[gid] - used_charge_sum[gid]);
                        
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        
                        pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                        //pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                    }
                    curcharge -= tmp as usize;
                }
                else {
                    if curcharge <= 24600 {
                        tmp = std::cmp::min((st.ev.c_ev_max - curcharge) as i32, -need[gid] - used_charge_sum[gid]);
                        
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                        //pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                    }
                    curcharge += tmp as usize;
                }
                if tmp != 0 {
                    pre = ev[j];
                }
                else {
                    cur_step -= st.graph.costs[pre][ev[j]];
                    cost -= nxcost;
                    curcharge += nxcost;
                }
            }
        }
        bestscore += cost as f64;
    }
    for gid in 0..st.nano_grid.n_grid {
        let mut add_pena = 0;
        let mut pena_sum = 0;
        if need[gid] >= 0 {
            for i in 0..lim_t {
                pena_sum += pena_list[gid][i];
                add_pena += if -p_pena[gid][i] > pena_sum {-p_pena[gid][i] - pena_sum} else {0}; 
            }
            grid_penalty[gid] = add_pena as f64 * 2.;
            grid_penalty[gid] += (need[gid] - used_charge_sum[gid]) as f64 / need_weight;
        }
        else {
            for i in 0..lim_t {
                pena_sum += pena_list[gid][i];
                add_pena += if p_pena[gid][i] > pena_sum {p_pena[gid][i] - pena_sum} else {0}; 
            }
            grid_penalty[gid] = add_pena as f64;
            grid_penalty[gid] -= used_charge_sum[gid] as f64 / avail_weight;
        }
        bestscore += grid_penalty[gid];
    }
    let mut cand_grid = Vec::new();
    for gid in 0..st.nano_grid.n_grid {
        let cnt = cntgrid[gid];
        let mut c;
        if need[gid] > 0 {
            c = (need[gid] + 12500 - 1) / 12500;
        }
        else {
            c = (-need[gid] + 15000 - 1) / 15000;
        }
        c = c.min(2);
        for _ in cnt..c {cand_grid.push(st.nano_grid.x[gid]);}
    }
    for _ in 0..3 {cand_grid.push(V);}
    let mut curscore = bestscore;
    let mut best = perm.clone();
    let starttemp = 250.0;
    let endtemp = 0.001;
    let num_itr = 200000;
    let invtl = 1.0 / num_itr as f64;

    if OUT_INFO {eprintln!("initscore: {}", curscore);}
    let selr = if cand_grid.len() != 0 {5} else {4};
    for lp in 0..num_itr {

        let ts = starttemp + (endtemp - starttemp) * lp as f64 * invtl;

        for eid1 in 0..st.ev.n_ev {
            let sel = st.rnd.next_int()%selr;
            let (mut id1, mut id2);
            let mut eid2 = 0;
            if sel <= 1 {
                id1 = st.rnd.nextn(5) + 1;
                id2 = st.rnd.nextn(5) + 1;
                if id1 == id2 {id2 = id2 % 5 + 2;}
                if id1 > id2 {swap(&mut id1, &mut id2);}
                if perm[eid1][id1] == V && perm[eid1][id2] == V {continue;}
                perm[eid1] = perm[eid1][..id1].into_iter().chain(perm[eid1][id1..id2+1].into_iter().rev()).chain(perm[eid1][id2+1..].into_iter()).copied().collect();
            }
            else if sel <= 3 {
                eid2 = st.rnd.nextn(st.ev.n_ev);
                while eid1 == eid2 {
                    eid2 = st.rnd.nextn(st.ev.n_ev);
                }
                id1 = st.rnd.nextn(5) + 1;
                id2 = st.rnd.nextn(5) + 1;
                if perm[eid1][id1] == V && perm[eid2][id2] == V {continue;}
                let tmp = perm[eid1][id1];
                perm[eid1][id1] = perm[eid2][id2];
                perm[eid2][id2] = tmp;
            }
            else {
                id1 = st.rnd.nextn(5) + 1;
                id2 = st.rnd.nextn(cand_grid.len());
                let tmp = perm[eid1][id1];
                perm[eid1][id1] = cand_grid[id2];
                cand_grid[id2] = tmp;
            }
            let mut nxscore = 0.;
            for i in used_charge_sum.iter_mut() {*i=0;}
            for i in pena_list.iter_mut().flat_map(|x| x) {*i=0;}
            for (eid, ev) in perm.iter().enumerate() {
                let mut curcharge = st.ev.charge[eid];
                let mut cost = st.ev.next[eid].1 * st.ev.delta_ev_move;
                curcharge -= cost;
                let mut pre = st.ev.next[eid].0;
                let mut cur_step = st.ev.next[eid].1;
                for j in 1..ev.len() {
                    let nxcost = st.graph.costs[pre][ev[j]] * st.ev.delta_ev_move;
                    if curcharge < nxcost {break;}
                    //cost += nxcost;
                    //curcharge -= nxcost;
                    if ev[j] != V {
                        cost += nxcost;
                        curcharge -= nxcost;
                        cur_step += st.graph.costs[pre][ev[j]];
                        let mut tmp = 0;
                        let gid = st.nano_grid.pos_to_id[ev[j]];
                        if need[gid] >= 0 {
                            if curcharge >= 1000 {
                                tmp = std::cmp::min(curcharge as i32 - rem, need[gid] - used_charge_sum[gid]);
                                
                                if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                                    tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                                }
                                else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                                pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                                used_charge_sum[gid] += tmp;
                                cur_step += tmp as usize / st.ev.v_ev_max;
                                //pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                            }
                            curcharge -= tmp as usize;
                        }
                        else {
                            if curcharge <= 24600 {
                                tmp = std::cmp::min((st.ev.c_ev_max - curcharge) as i32, -need[gid] - used_charge_sum[gid]);
                                
                                if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                                    tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                                }
                                else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                                pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                                used_charge_sum[gid] += tmp;
                                cur_step += tmp as usize / st.ev.v_ev_max;
                                //pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                            }
                            curcharge += tmp as usize;
                        }
                        if tmp != 0 {
                            pre = ev[j];
                        }
                        else {
                            cur_step -= st.graph.costs[pre][ev[j]];
                            cost -= nxcost;
                            curcharge += nxcost;
                        }
                    }
                }
                nxscore += cost as f64;
            }
            for gid in 0..st.nano_grid.n_grid {
                let mut add_pena = 0;
                let mut pena_sum = 0;
                if need[gid] >= 0 {
                    for i in 0..lim_t {
                        pena_sum += pena_list[gid][i];
                        add_pena += if -p_pena[gid][i] > pena_sum {-p_pena[gid][i] - pena_sum} else {0}; 
                    }
                    grid_penalty[gid] = add_pena as f64 * 2.;
                    grid_penalty[gid] += (need[gid] - used_charge_sum[gid]) as f64 / need_weight;
                }
                else {
                    for i in 0..lim_t {
                        pena_sum += pena_list[gid][i];
                        add_pena += if p_pena[gid][i] > pena_sum {p_pena[gid][i] - pena_sum} else {0}; 
                    }
                    grid_penalty[gid] = add_pena as f64;
                    grid_penalty[gid] -= used_charge_sum[gid] as f64 / avail_weight;
                }
                nxscore += grid_penalty[gid];
            }

            let diff = curscore - nxscore;
            if diff >= 0. || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                curscore -= diff;
                if bestscore > curscore {
                    bestscore = curscore;
                    best = perm.clone();
                }
            }
            else {
                if sel <= 1 {
                    perm[eid1] = perm[eid1][..id1].into_iter().chain(perm[eid1][id1..id2+1].into_iter().rev()).chain(perm[eid1][id2+1..].into_iter()).copied().collect();
                }
                else if sel <= 3 {
                    let tmp = perm[eid1][id1];
                    perm[eid1][id1] = perm[eid2][id2];
                    perm[eid2][id2] = tmp;
                }
                else {
                    let tmp = perm[eid1][id1];
                    perm[eid1][id1] = cand_grid[id2];
                    cand_grid[id2] = tmp;
                }
            }
        }
    }
    if OUT_INFO {eprintln!("afterscore: {}", bestscore);}
    //perm bestから再構築
    for i in used_charge_sum.iter_mut() {*i=0;}
    for i in pena_list.iter_mut().flat_map(|x| x) {*i=0;}
    let mut nxscore = 0.;
    for (eid, ev) in best.iter().enumerate() {
        perm[eid].clear();
        let mut curcharge = st.ev.charge[eid];
        let mut cost = st.ev.next[eid].1 * st.ev.delta_ev_move;
        curcharge -= cost;
        let mut pre = st.ev.next[eid].0;
        perm[eid].push(pre);
        let mut cur_step = st.ev.next[eid].1;
        for j in 1..ev.len() {
            let nxcost = st.graph.costs[pre][ev[j]] * st.ev.delta_ev_move;
            if curcharge < nxcost {break;}
            //cost += nxcost;
            //curcharge -= nxcost;
            if ev[j] != V {
                cost += nxcost;
                curcharge -= nxcost;
                cur_step += st.graph.costs[pre][ev[j]];
                let mut tmp = 0;
                let gid = st.nano_grid.pos_to_id[ev[j]];
                if need[gid] >= 0 {
                    if curcharge >= 1000 {
                        tmp = std::cmp::min(curcharge as i32 - rem, need[gid] - used_charge_sum[gid]);
                        
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                        //pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                    }
                    curcharge -= tmp as usize;
                }
                else {
                    if curcharge < 24600 {
                        tmp = std::cmp::min((st.ev.c_ev_max - curcharge) as i32, -need[gid] - used_charge_sum[gid]);
                        
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                        //pena_list[gid][3.min(cur_step/DIV_RANGE)] += tmp;
                    }
                    curcharge += tmp as usize;
                }
                if tmp != 0 {
                    pre = ev[j];
                    perm[eid].push(pre);
                }
                else {
                    cur_step -= st.graph.costs[pre][ev[j]];
                    cost -= nxcost;
                    curcharge += nxcost;
                }
            }
        }
        nxscore += cost as f64;
        perm[eid].push(V);
    }
    for gid in 0..st.nano_grid.n_grid {
        let mut add_pena = 0;
        let mut pena_sum = 0;
        if need[gid] >= 0 {
            for i in 0..lim_t {
                pena_sum += pena_list[gid][i];
                add_pena += if -p_pena[gid][i] > pena_sum {-p_pena[gid][i] - pena_sum} else {0}; 
            }
            grid_penalty[gid] = add_pena as f64 * 2.;
            grid_penalty[gid] += (need[gid] - used_charge_sum[gid]) as f64 / need_weight;
        }
        else {
            for i in 0..lim_t {
                pena_sum += pena_list[gid][i];
                add_pena += if p_pena[gid][i] > pena_sum {p_pena[gid][i] - pena_sum} else {0}; 
            }
            grid_penalty[gid] = add_pena as f64;
            grid_penalty[gid] -= used_charge_sum[gid] as f64 / avail_weight;
        }
        nxscore += grid_penalty[gid];
    }
    if OUT_INFO {eprintln!("afterscore: {} / {}", bestscore, nxscore);}
    perm
}

fn optimize_ev_route2 (st: &mut State, mut perm: Vec<Vec<usize>>, need: &Vec<i32>, div_t: usize) -> Vec<Vec<usize>> {
    let mut bestscore = 0.;
    let mut grid_penalty = vec![0.; st.nano_grid.n_grid];
    let mut used_charge_sum = vec![0; st.nano_grid.n_grid];
    let (need_weight, avail_weight) = (3., 10.);
    let rem = if div_t <= 12 {1000} else {0};
    let lim_step = match st.nano_grid.day_type {
        DayType::FineSquall => {50},
        _ => {300}//inf
    };
    for (eid, ev) in perm.iter().enumerate() {
        let mut curcharge = st.ev.charge[eid];
        let mut cost = st.ev.next[eid].1 * st.ev.delta_ev_move;
        curcharge -= cost;
        let mut pre = st.ev.next[eid].0;
        let mut cur_step = st.ev.next[eid].1;
        for j in 1..ev.len() {
            let nxcost = st.graph.costs[pre][ev[j]] * st.ev.delta_ev_move;
            if curcharge < nxcost {break;}
            //cost += nxcost;
            //curcharge -= nxcost;
            if ev[j] != V {
                cost += nxcost;
                curcharge -= nxcost;
                cur_step += st.graph.costs[pre][ev[j]];
                let mut tmp = 0;
                let gid = st.nano_grid.pos_to_id[ev[j]];
                if need[gid] >= 0 {
                    if curcharge >= 1000 {
                        tmp = std::cmp::min(curcharge as i32 - rem, need[gid] - used_charge_sum[gid]);
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                    }
                    curcharge -= tmp as usize;
                }
                else {
                    if curcharge <= 24600 {
                        tmp = std::cmp::min((st.ev.c_ev_max - curcharge) as i32, -need[gid] - used_charge_sum[gid]);
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                    }
                    curcharge += tmp as usize;
                }
                if tmp != 0 {
                    pre = ev[j];
                }
                else {
                    cur_step -= st.graph.costs[pre][ev[j]];
                    cost -= nxcost;
                    curcharge += nxcost;
                }
            }
        }
        bestscore += cost as f64;
    }
    for gid in 0..st.nano_grid.n_grid {
        if need[gid] >= 0 {
            let pr = if st.nano_grid.in_sudden[gid] {2.} else {1.};
            grid_penalty[gid] += (need[gid] - used_charge_sum[gid]) as f64 / need_weight * pr;
        }
        else {
            grid_penalty[gid] -= used_charge_sum[gid] as f64 / avail_weight;
        }
        bestscore += grid_penalty[gid];
    }
    let mut curscore = bestscore;
    let mut best = perm.clone();
    let starttemp = 250.0;
    let endtemp = 0.001;
    let num_itr = 200000;
    let invtl = 1.0 / num_itr as f64;

    if OUT_INFO {eprintln!("initscore: {}", curscore);}
    for lp in 0..num_itr {

        let ts = starttemp + (endtemp - starttemp) * lp as f64 * invtl;

        for eid1 in 0..st.ev.n_ev {
            let sel = st.rnd.next_int()&1;
            let (mut id1, mut id2);
            let mut eid2 = 0;
            if sel == 0 {
                id1 = st.rnd.nextn(5) + 1;
                id2 = st.rnd.nextn(5) + 1;
                if id1 == id2 {id2 = id2 % 5 + 2;}
                if id1 > id2 {swap(&mut id1, &mut id2);}
                if perm[eid1][id1] == V && perm[eid1][id2] == V {continue;}
                perm[eid1] = perm[eid1][..id1].into_iter().chain(perm[eid1][id1..id2+1].into_iter().rev()).chain(perm[eid1][id2+1..].into_iter()).copied().collect();
            }
            else {
                eid2 = st.rnd.nextn(st.ev.n_ev);
                while eid1 == eid2 {
                    eid2 = st.rnd.nextn(st.ev.n_ev);
                }
                id1 = st.rnd.nextn(5) + 1;
                id2 = st.rnd.nextn(5) + 1;
                if perm[eid1][id1] == V && perm[eid2][id2] == V {continue;}
                let tmp = perm[eid1][id1];
                perm[eid1][id1] = perm[eid2][id2];
                perm[eid2][id2] = tmp;
            }
            let mut nxscore = 0.;
            for i in used_charge_sum.iter_mut() {*i=0;}
            for i in grid_penalty.iter_mut() {*i=0.;}

            for (eid, ev) in perm.iter().enumerate() {
                let mut curcharge = st.ev.charge[eid];
                let mut cost = st.ev.next[eid].1 * st.ev.delta_ev_move;
                curcharge -= cost;
                let mut pre = st.ev.next[eid].0;
                let mut cur_step = st.ev.next[eid].1;
                for j in 1..ev.len() {
                    let nxcost = st.graph.costs[pre][ev[j]] * st.ev.delta_ev_move;
                    if curcharge < nxcost {break;}
                    //cost += nxcost;
                    //curcharge -= nxcost;
                    if ev[j] != V {
                        cost += nxcost;
                        curcharge -= nxcost;
                        cur_step += st.graph.costs[pre][ev[j]];
                        let mut tmp = 0;
                        let gid = st.nano_grid.pos_to_id[ev[j]];
                        if need[gid] >= 0 {
                            if curcharge >= 1000 {
                                tmp = std::cmp::min(curcharge as i32 - rem, need[gid] - used_charge_sum[gid]);
                                if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                                    tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                                }
                                else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                                used_charge_sum[gid] += tmp;
                                cur_step += tmp as usize / st.ev.v_ev_max;
                            }
                            curcharge -= tmp as usize;
                        }
                        else {
                            if curcharge <= 24600 {
                                tmp = std::cmp::min((st.ev.c_ev_max - curcharge) as i32, -need[gid] - used_charge_sum[gid]);
                                if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                                    tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                                }
                                else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                                used_charge_sum[gid] += tmp;
                                cur_step += tmp as usize / st.ev.v_ev_max;
                            }
                            curcharge += tmp as usize;
                        }
                        if tmp != 0 {
                            pre = ev[j];
                        }
                        else {
                            cur_step -= st.graph.costs[pre][ev[j]];
                            cost -= nxcost;
                            curcharge += nxcost;
                        }
                    }
                }
                nxscore += cost as f64;
            }
            for gid in 0..st.nano_grid.n_grid {
                if need[gid] >= 0 {
                    let pr = if st.nano_grid.in_sudden[gid] {2.} else {1.};
                    grid_penalty[gid] += (need[gid] - used_charge_sum[gid]) as f64 / need_weight * pr;
                }
                else {
                    grid_penalty[gid] -= used_charge_sum[gid] as f64 / avail_weight;
                }
                nxscore += grid_penalty[gid];
            }

            let diff = curscore - nxscore;
            if diff >= 0. || diff as f64 > st.ln[st.rnd.next_int() & 65535] * ts {
                curscore -= diff;
                if bestscore > curscore {
                    bestscore = curscore;
                    best = perm.clone();
                }
            }
            else {
                if sel == 0 {
                    perm[eid1] = perm[eid1][..id1].into_iter().chain(perm[eid1][id1..id2+1].into_iter().rev()).chain(perm[eid1][id2+1..].into_iter()).copied().collect();
                }
                else {
                    let tmp = perm[eid1][id1];
                    perm[eid1][id1] = perm[eid2][id2];
                    perm[eid2][id2] = tmp;
                }
            }
        }
    }
    if OUT_INFO {eprintln!("afterscore: {}", bestscore);}
    //perm bestから再構築
    for i in used_charge_sum.iter_mut() {*i=0;}
    for i in grid_penalty.iter_mut() {*i=0.;}
    let mut nxscore = 0.;
    for (eid, ev) in best.iter().enumerate() {
        perm[eid].clear();
        let mut curcharge = st.ev.charge[eid];
        let mut cost = st.ev.next[eid].1 * st.ev.delta_ev_move;
        curcharge -= cost;
        let mut pre = st.ev.next[eid].0;
        perm[eid].push(pre);
        let mut cur_step = st.ev.next[eid].1;
        for j in 1..ev.len() {
            let nxcost = st.graph.costs[pre][ev[j]] * st.ev.delta_ev_move;
            if curcharge < nxcost {break;}
            //cost += nxcost;
            //curcharge -= nxcost;
            if ev[j] != V {
                cost += nxcost;
                curcharge -= nxcost;
                cur_step += st.graph.costs[pre][ev[j]];
                let mut tmp = 0;
                let gid = st.nano_grid.pos_to_id[ev[j]];
                if need[gid] >= 0 {
                    if curcharge >= 1000 {
                        tmp = std::cmp::min(curcharge as i32 - rem, need[gid] - used_charge_sum[gid]);
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                    }
                    curcharge -= tmp as usize;
                }
                else {
                    if curcharge < 24600 {
                        tmp = std::cmp::min((st.ev.c_ev_max - curcharge) as i32, -need[gid] - used_charge_sum[gid]);
                        if lim_step > cur_step && st.nano_grid.in_sudden[gid] {
                            tmp = tmp.min(((lim_step - cur_step) * st.ev.v_ev_max) as i32);
                        }
                        else if st.nano_grid.in_sudden[gid] {tmp = 0;}
                        used_charge_sum[gid] += tmp;
                        cur_step += tmp as usize / st.ev.v_ev_max;
                    }
                    curcharge += tmp as usize;
                }
                if tmp != 0 {
                    pre = ev[j];
                    perm[eid].push(pre);
                }
                else {
                    cur_step -= st.graph.costs[pre][ev[j]];
                    cost -= nxcost;
                    curcharge += nxcost;
                }
            }
        }
        nxscore += cost as f64;
        perm[eid].push(V);
    }
    for gid in 0..st.nano_grid.n_grid {
        if need[gid] >= 0 {
            let pr = if st.nano_grid.in_sudden[gid] {2.} else {1.};
            grid_penalty[gid] += (need[gid] - used_charge_sum[gid]) as f64 / need_weight * pr;
        }
        else {
            grid_penalty[gid] -= used_charge_sum[gid] as f64 / avail_weight;
        }
        nxscore += grid_penalty[gid];
    }
    if OUT_INFO {eprintln!("afterscore: {} / {}", bestscore, nxscore);}
    perm
}

fn construct_ev_route (st: &mut State, need: Vec<(i32, usize)>, div_t: usize) -> (Vec<Vec<usize>>, Vec<i32>) {
    let mut perm: Vec<Vec<usize>> = (0..st.ev.n_ev).map(|i| vec![st.ev.next[i].0]).collect();
    let mut used_grid: Vec<usize> = vec![0; st.nano_grid.n_grid];
    
    for i in 0..st.nano_grid.n_grid {
        st.nano_grid.visit_ev_id[i].clear();
        st.nano_grid.chrg_end_time[i].clear();
        let f_t = st.nano_grid.n_div.min(div_t + 1); 
        let pat = st.nano_grid.patterns[i];
        let mut overflow = st.nano_grid.charge[i] as i32 + st.nano_grid.pw_expected[pat][f_t] - st.nano_grid.pw_expected[pat][div_t];
        if st.nano_grid.in_sudden[i] {
            overflow += match st.nano_grid.day_type {
                DayType::FineSquall => {-50000},
                DayType::RainyFine => {50000},
                _ => {0}
            }
        }
        st.nano_grid.overflow[i] = overflow - st.nano_grid.c_grid_max as i32;
        if st.nano_grid.overflow[i] > 0 {st.complete = false;}
    }
    let flag = match st.nano_grid.day_type {
        DayType::FineSquall | DayType::RainyFine => {true},
        _ => {false}
    };
    
    if st.sudden && flag {
        for i in 0..st.nano_grid.n_grid {
            let j = st.nano_grid.pos_to_id[need[i].1];
            let pat = st.nano_grid.patterns[j];
            let tmp = st.nano_grid.charge[j] as i32 + st.nano_grid.pw_expected[pat][20.min(div_t+3)] - st.nano_grid.pw_expected[pat][div_t];
            if st.nano_grid.in_sudden[j] && need[i].0 > 0 {
                used_grid[i] = 2;
            }
            else if tmp < 0 || st.nano_grid.charge[j] == 0 {
                used_grid[i] = 1;
            }
            else {used_grid[i] = 0;}
        }
    }
    else {
        for i in 0..st.nano_grid.n_grid {
            let j = st.nano_grid.pos_to_id[need[i].1];
            if need[i].0 <= 500 && st.nano_grid.charge[j] != 0 {continue;}
            let cnt = std::cmp::max(0, ((need[i].0 + 12500 - 1) / 12500) as usize);
            used_grid[i] = std::cmp::min(2, cnt); 
        }
    }
    let co = match st.nano_grid.day_type {
        DayType::Fine | DayType::FineSquall => {5000},
        DayType::Rainy | DayType::RainyFine => {1000}
    };
    let mut gid = 0;
    while used_grid[gid] != 0 {
        let mut mn = INF;
        let mut eid = 0;
        let nid = need[gid].1;
        for ev in 0..st.ev.n_ev {
            let plen = perm[ev].len();
            let chr_amnt = st.ev.charge[ev];
            
            if plen > 3 || plen > chr_amnt / co {continue;}
            let &posid = perm[ev].last().unwrap();
            if (plen == 1 || plen > 1 && posid != nid) && mn > st.graph.costs[posid][nid] {
                eid = ev;
                mn = st.graph.costs[posid][nid];
            }
        }
        if mn == INF {break;}
        used_grid[gid] -= 1;
        perm[eid].push(nid);
        st.nano_grid.visit_ev_id[gid].push_back(eid);
        while gid < need.len() && used_grid[gid] == 0 {gid += 1;}
        if gid == need.len() {break;}
    }
    let mut mxpat = 12500;
    for pat in 0..st.nano_grid.n_pattern {
        mxpat = mxpat.max(st.nano_grid.pw_expected[pat][st.nano_grid.n_div] + st.nano_grid.c_grid_init as i32);
    }
    let limit = match st.nano_grid.day_type {
        DayType::RainyFine => {-mxpat + 5000},
        DayType::FineSquall => {-500},
        _ => {-3000}
    };

    let flag = match st.nano_grid.day_type {
        DayType::RainyFine => {true},
        _ => {false}
    };
    
    if st.sudden && flag {
        for i in 0..st.nano_grid.n_grid {
            let j = st.nano_grid.pos_to_id[need[i].1];
            let pat = st.nano_grid.patterns[j];
            let tmp = st.nano_grid.charge[j] as i32 + st.nano_grid.pw_expected[pat][20.min(div_t+3)] - st.nano_grid.pw_expected[pat][div_t];
            if tmp > st.nano_grid.c_grid_max as i32 || st.nano_grid.in_sudden[j] && need[i].0 < 0 {
                used_grid[i] = 2;
            }
            else {used_grid[i] = 0;}
        }
    }
    else {
        for i in 0..st.nano_grid.n_grid {
            let j = st.nano_grid.pos_to_id[need[i].1];
            if need[i].0 >= limit || st.nano_grid.charge[j] < 500 {
                used_grid[i] = 0;
                continue;
            }
            let cnt = std::cmp::max(0, ((-need[i].0 + 15000 - 1) / 15000) as usize);//hyper_param
            used_grid[i] = std::cmp::min(2, cnt); 
            if st.nano_grid.overflow[j] > 0 {
                used_grid[i] = 2;
            }
        }
    }
    let border = 22500;
    let mut gid = need.len() - 1;
    while gid > 0 && st.nano_grid.visit_ev_id[gid].len() == 2 {gid -= 1;}
    while gid > 0 && used_grid[gid] != 0 {
        let mut mn = INF;
        let mut eid = 0;
        let nid = need[gid].1;
        let j = st.nano_grid.pos_to_id[need[gid].1];
        for ev in 0..st.ev.n_ev {
            if st.ev.charge[ev] >= border || st.complete && st.nano_grid.overflow[j] < 500 {continue;}
            if perm[ev].len() > 5 {continue;}
            let &posid = perm[ev].last().unwrap();
            if  mn > st.graph.costs[posid][nid] && st.ev.charge[ev] >= st.graph.costs[posid][nid] * st.ev.delta_ev_move {
                eid = ev;
                mn = st.graph.costs[posid][nid];
            }
        }
        if mn == INF {break;}
        used_grid[gid] -= 1;
        perm[eid].push(nid);
        st.nano_grid.visit_ev_id[gid].push_back(eid);
        if used_grid[gid] == 0 {gid -= 1;}
        while gid > 0 && st.nano_grid.visit_ev_id[gid].len() == 2 {gid -= 1;}
        if gid == 0 {break;}
    }
    for i in 0..st.ev.n_ev {
        while perm[i].len() != 7 {perm[i].push(V);}
    }

    let mut rneed = vec![0i32; st.nano_grid.n_grid];
    for &(nd, i) in need.iter() {
        rneed[st.nano_grid.pos_to_id[i]] = nd;
    }

    //let perm = optimize_ev_route(st, perm, &rneed, div_t);
    
    let perm = match st.nano_grid.day_type {
        DayType::Fine | DayType::Rainy => {optimize_ev_route(st, perm, &rneed, div_t)},
        DayType::FineSquall | DayType::RainyFine =>  {
            if st.sudden {optimize_ev_route2(st, perm, &rneed, div_t)}
            else {optimize_ev_route(st, perm, &rneed, div_t)}
        }
    };
    
    
    (perm, rneed)
}

fn construct_commands(st: &mut State, perm: &Vec<Vec<usize>>, need: &mut Vec<i32>, tt: usize, div_t: usize) {
    for i in st.nano_grid.use_count.iter_mut() {*i = 0;}
    
    for i in 0..st.ev.n_ev {
        let mut curcharge = st.ev.charge[i];
        let mut pos = st.ev.next[i].0; 
        st.ev.commands[i].clear();
        for _ in 0..st.ev.next[i].1 {
            st.ev.commands[i].push_back(Command::Move(pos));
            curcharge -= st.ev.delta_ev_move;
        }
        for &nx in perm[i].iter().skip(1) {
            if nx == V {break;}
            let nxid = st.nano_grid.pos_to_id[nx];
            if st.graph.costs[pos][nx] * st.ev.delta_ev_move > curcharge {continue;}

            //if need[nxid].abs() < 1500 {continue;}
            while pos != nx {
                let subnx = st.graph.next[pos][nx];
                let dist = st.graph.costs[pos][subnx];
                for _ in 0..dist {
                    st.ev.commands[i].push_back(Command::Move(subnx));
                    curcharge -= st.ev.delta_ev_move;
                }
                pos = subnx;
            }
            let nxid = st.nano_grid.pos_to_id[nx];
            
            st.nano_grid.use_count[nxid] += 1;
            if need[nxid] < 0 {
                //if curcharge >= 24600 {break;}
                let charge_turn = std::cmp::min(st.ev.c_ev_max - curcharge, (-need[nxid]) as usize) / st.ev.v_ev_max;
                for _ in 0..charge_turn {
                    let vc = std::cmp::min((-need[nxid]) as usize, st.ev.v_ev_max);
                    st.ev.commands[i].push_back(Command::ChargeFromGrid(vc));
                    need[nxid] += vc as i32;
                    curcharge += vc;
                    if curcharge >= 24600 {break;}
                }
                if charge_turn > 0 {st.nano_grid.chrg_end_time[nxid].push_back(tt + st.ev.commands[i].len());}
            }
            else {
                let rem = if div_t <= 12 {1000} else {0};
                if curcharge <= rem {continue;}
                let charge_turn = std::cmp::min(curcharge - rem, need[nxid] as usize + st.ev.v_ev_max - 1) / st.ev.v_ev_max;
                for _ in 0..charge_turn {
                    let vc = std::cmp::min(need[nxid] as usize, st.ev.v_ev_max);
                    st.ev.commands[i].push_back(Command::ChargeToGrid(vc));
                    need[nxid] -= vc as i32;
                    curcharge -= vc;
                }
                if charge_turn > 0 {st.nano_grid.chrg_end_time[nxid].push_back(tt + st.ev.commands[i].len());}
            }
            //st.nano_grid.use_count[nxid] += 1;
        }
    }
}

fn construct_additional_commands(st: &mut State, need: &mut Vec<i32>, eid: usize, tt: usize, div_t: usize) {
    let mut nxid = V;
    let mut pos = st.ev.next[eid].0;
    let acost = match st.ev.pre_commands[eid] {
        Command::Stay => {st.ev.next[eid].1},
        Command::Move(_to) => {
            if st.ev.next[eid].1 == 0 {1} else {st.ev.next[eid].1 - 1}
        },
        _ => {0}
    };
    
    let mut curcharge = {
        match st.ev.pre_commands[eid] {
            Command::Stay => {st.ev.charge[eid]},
            Command::Move(_to) => {
                st.ev.charge[eid] - st.ev.delta_ev_move
            },
            Command::ChargeFromGrid(c) => {
                std::cmp::min(st.ev.c_ev_max, st.ev.charge[eid] + c)
            },
            Command::ChargeToGrid(c) => {
                st.ev.charge[eid] - c
            },
        }
    };

    let init_dist = match st.nano_grid.day_type {
        DayType::Fine => {3},
        DayType::RainyFine => {0},
        DayType::Rainy => {3},
        DayType::FineSquall => {6}
    };
    
    if curcharge >= 1000 {
        let mut mn = init_dist;
        for i in 0..st.nano_grid.n_grid {
            let cost = st.graph.costs[pos][st.nano_grid.x[i]] + acost;
            if st.nano_grid.chrg_end_time[i].len() == 2 {
                let &fr = st.nano_grid.chrg_end_time[i].front().unwrap();
                let &bk = st.nano_grid.chrg_end_time[i].back().unwrap();
                if tt + cost < std::cmp::min(fr, bk) {continue;}
            }
            if (st.nano_grid.charge[i] == 0 || need[i] > 0 && st.nano_grid.in_sudden[i])
            && cost * st.ev.delta_ev_move + 400 <= curcharge && mn > cost {
                mn = cost;
                nxid = i;
            }
        }
    }
    
    let init_dist = match st.nano_grid.day_type {
        DayType::Fine => {10},
        DayType::RainyFine => {6},
        DayType::Rainy => {6},
        DayType::FineSquall => {10}
    };
    let bd = if st.sudden {
        match st.nano_grid.day_type {
            DayType::RainyFine => {22000},
            _ => {22000}
        }
    } else {22000}; 
    if nxid == V && curcharge < bd {
        let mut mn = init_dist;
        for i in 0..st.nano_grid.n_grid {
            let cost = st.graph.costs[pos][st.nano_grid.x[i]] + acost;
            if st.nano_grid.chrg_end_time[i].len() == 2 {
                let &fr = st.nano_grid.chrg_end_time[i].front().unwrap();
                let &bk = st.nano_grid.chrg_end_time[i].back().unwrap();
                if tt + cost < std::cmp::min(fr, bk) {continue;}
            }
            if (st.nano_grid.charge[i] == st.nano_grid.c_grid_max || st.nano_grid.overflow[i] > 1000)
            && cost * st.ev.delta_ev_move <= curcharge && mn > cost {
                mn = cost;
                nxid = i;
            }
        }
    }
    
    let init_dist = match st.nano_grid.day_type {
        DayType::Fine => {0},
        DayType::RainyFine => {3},
        DayType::Rainy => {3},
        DayType::FineSquall => {0}
    };
    if nxid == V && curcharge >= 2000 {
        let mut mn = init_dist;
        for i in 0..st.nano_grid.n_grid {
            let cost = st.graph.costs[pos][st.nano_grid.x[i]] + acost;
            if st.nano_grid.chrg_end_time[i].len() == 2 {
                let &fr = st.nano_grid.chrg_end_time[i].front().unwrap();
                let &bk = st.nano_grid.chrg_end_time[i].back().unwrap();
                if tt + cost < std::cmp::min(fr, bk) {continue;}
            }
            if need[i] >= 1000 && cost * st.ev.delta_ev_move + 400 <= curcharge && mn > cost {
                mn = cost;
                nxid = i;
            }
        }
    }
    
    let last_count = match st.nano_grid.day_type {
        DayType::FineSquall => {1},
        _ => {8}
    };
    
    let init_dist = match st.nano_grid.day_type {
        DayType::Fine => {6},
        DayType::RainyFine => {6},
        DayType::Rainy => {6},
        DayType::FineSquall => {6}
    };
    if nxid == V && !st.complete && curcharge >= 500 && (div_t + last_count >= st.nano_grid.n_div) {
        let mut mn = init_dist;
        for i in 0..st.nano_grid.n_grid {
            let cost = st.graph.costs[pos][st.nano_grid.x[i]] + acost;
            if st.nano_grid.chrg_end_time[i].len() == 2 {
                let &fr = st.nano_grid.chrg_end_time[i].front().unwrap();
                let &bk = st.nano_grid.chrg_end_time[i].back().unwrap();
                if tt + cost < std::cmp::min(fr, bk) {continue;}
            }
            if need[i] >= 400 && cost * st.ev.delta_ev_move + 500 <= curcharge && mn > cost {
                mn = cost;
                nxid = i;
            }
        }
    }
    /*
    let init_dist = match st.nano_grid.day_type {
        DayType::Fine => {0},
        DayType::RainyFine => {0},
        DayType::Rainy => {0},
        DayType::FineSquall => {20}
    };
    if nxid == V && !st.complete && curcharge >= 400 && (div_t + last_count >= st.nano_grid.n_div) && (div_t + 1 < st.nano_grid.n_div) {
        let mut mn = init_dist;
        for i in 0..st.nano_grid.n_grid {
            let cost = st.graph.costs[pos][st.nano_grid.x[i]] + acost;
            if st.nano_grid.chrg_end_time[i].len() == 2 {
                let &fr = st.nano_grid.chrg_end_time[i].front().unwrap();
                let &bk = st.nano_grid.chrg_end_time[i].back().unwrap();
                if tt + cost < std::cmp::min(fr, bk) {continue;}
            }
            if need[i] >= 1000 && st.nano_grid.charge[i] == 0 && cost * st.ev.delta_ev_move + 400 <= curcharge && mn > cost {
                mn = cost;
                nxid = i;
            }
        }
    }
    */
    if nxid != V {
        st.ev.commands[eid].clear();//should be empty
        let gid = st.nano_grid.x[nxid];
        for _ in 0..acost {
            st.ev.commands[eid].push_back(Command::Move(pos));
            curcharge -= st.ev.delta_ev_move;
        }
        while pos != gid {
            let subnx = st.graph.next[pos][gid];
            let dist = st.graph.costs[pos][subnx];
            for _ in 0..dist {
                st.ev.commands[eid].push_back(Command::Move(subnx));
                curcharge -= st.ev.delta_ev_move;
            }
            pos = subnx;
        }
        if need[nxid] < 0 {
            let charge_turn = std::cmp::min(st.ev.c_ev_max - curcharge, (-need[nxid]) as usize) / st.ev.v_ev_max;
            for _ in 0..charge_turn {
                let vc = std::cmp::min((-need[nxid]) as usize, st.ev.v_ev_max);
                st.ev.commands[eid].push_back(Command::ChargeFromGrid(vc));
                need[nxid] += vc as i32;
                curcharge += vc;
                if curcharge >= 24500 {break;}
            }
        }
        else {
            //let mut rem = if div_t + last_count >= st.nano_grid.n_div {0} else {1200};
            let mut rem = if div_t + last_count >= st.nano_grid.n_div || st.nano_grid.charge[nxid] == 0 {0} else {1000};
            rem = rem.min(curcharge);
            let charge_turn = (std::cmp::min(curcharge - rem, need[nxid] as usize) + st.ev.v_ev_max - 1) / st.ev.v_ev_max;
            for _ in 0..charge_turn {
                let mut vc = std::cmp::min(need[nxid] as usize, st.ev.v_ev_max);
                vc = vc.min(curcharge);
                st.ev.commands[eid].push_back(Command::ChargeToGrid(vc));
                need[nxid] -= vc as i32;
                curcharge -= vc;
            }
        }

        if st.nano_grid.chrg_end_time[nxid].len() == 2 { 
            let &fr = st.nano_grid.chrg_end_time[nxid].front().unwrap();
            let &bk = st.nano_grid.chrg_end_time[nxid].back().unwrap();
            if fr < bk {
                st.nano_grid.chrg_end_time[nxid].pop_front();
            }
            else {
                st.nano_grid.chrg_end_time[nxid].pop_back();
            }
        }
        st.nano_grid.chrg_end_time[nxid].push_back(tt + st.ev.commands[eid].len());
        st.nano_grid.visit_ev_id[nxid].push_back(eid);
    }
}

fn construct_additional_commands2(st: &mut State, need: &mut Vec<i32>, eid: usize, tt: usize, div_t: usize) {
    let mut nxid = V;
    let mut pos = st.ev.next[eid].0;
    let acost = match st.ev.pre_commands[eid] {
        Command::Stay => {st.ev.next[eid].1},
        Command::Move(_to) => {
            if st.ev.next[eid].1 == 0 {1} else {st.ev.next[eid].1 - 1}
        },
        _ => {0}
    };
    
    let mut curcharge = {
        match st.ev.pre_commands[eid] {
            Command::Stay => {st.ev.charge[eid]},
            Command::Move(_to) => {
                st.ev.charge[eid] - st.ev.delta_ev_move
            },
            Command::ChargeFromGrid(c) => {
                std::cmp::min(st.ev.c_ev_max, st.ev.charge[eid] + c)
            },
            Command::ChargeToGrid(c) => {
                st.ev.charge[eid] - c
            },
        }
    };

    let last_count = match st.nano_grid.day_type {
        DayType::FineSquall => {2},
        _ => {8}
    };
    let init_dist = match st.nano_grid.day_type {
        DayType::Fine => {0},
        DayType::RainyFine => {0},
        DayType::Rainy => {0},
        DayType::FineSquall => {20}
    };
    if nxid == V && !st.complete && curcharge >= 400 && (div_t + last_count >= st.nano_grid.n_div) && (div_t + 1 < st.nano_grid.n_div) {
        let mut mn = init_dist;
        for i in 0..st.nano_grid.n_grid {
            let cost = st.graph.costs[pos][st.nano_grid.x[i]] + acost;
            if st.nano_grid.chrg_end_time[i].len() == 2 {
                let &fr = st.nano_grid.chrg_end_time[i].front().unwrap();
                let &bk = st.nano_grid.chrg_end_time[i].back().unwrap();
                if tt + cost < std::cmp::min(fr, bk) {continue;}
            }
            if need[i] >= 1000 && st.nano_grid.charge[i] == 0 && cost * st.ev.delta_ev_move + 400 <= curcharge && mn > cost {
                mn = cost;
                nxid = i;
            }
        }
    }

    if nxid != V {
        st.ev.commands[eid].clear();//should be empty
        let gid = st.nano_grid.x[nxid];
        for _ in 0..acost {
            st.ev.commands[eid].push_back(Command::Move(pos));
            curcharge -= st.ev.delta_ev_move;
        }
        while pos != gid {
            let subnx = st.graph.next[pos][gid];
            let dist = st.graph.costs[pos][subnx];
            for _ in 0..dist {
                st.ev.commands[eid].push_back(Command::Move(subnx));
                curcharge -= st.ev.delta_ev_move;
            }
            pos = subnx;
        }
        if need[nxid] < 0 {
            let charge_turn = std::cmp::min(st.ev.c_ev_max - curcharge, (-need[nxid]) as usize) / st.ev.v_ev_max;
            for _ in 0..charge_turn {
                let vc = std::cmp::min((-need[nxid]) as usize, st.ev.v_ev_max);
                st.ev.commands[eid].push_back(Command::ChargeFromGrid(vc));
                need[nxid] += vc as i32;
                curcharge += vc;
                if curcharge >= 24500 {break;}
            }
        }
        else {
            //let mut rem = if div_t + last_count >= st.nano_grid.n_div {0} else {1200};
            let mut rem = if div_t + last_count >= st.nano_grid.n_div || st.nano_grid.charge[nxid] == 0 {0} else {1000};
            rem = rem.min(curcharge);
            let charge_turn = (std::cmp::min(curcharge - rem, need[nxid] as usize) + st.ev.v_ev_max - 1) / st.ev.v_ev_max;
            for _ in 0..charge_turn {
                let mut vc = std::cmp::min(need[nxid] as usize, st.ev.v_ev_max);
                vc = vc.min(curcharge);
                st.ev.commands[eid].push_back(Command::ChargeToGrid(vc));
                need[nxid] -= vc as i32;
                curcharge -= vc;
            }
        }

        if st.nano_grid.chrg_end_time[nxid].len() == 2 { 
            let &fr = st.nano_grid.chrg_end_time[nxid].front().unwrap();
            let &bk = st.nano_grid.chrg_end_time[nxid].back().unwrap();
            if fr < bk {
                st.nano_grid.chrg_end_time[nxid].pop_front();
            }
            else {
                st.nano_grid.chrg_end_time[nxid].pop_back();
            }
        }
        st.nano_grid.chrg_end_time[nxid].push_back(tt + st.ev.commands[eid].len());
        st.nano_grid.visit_ev_id[nxid].push_back(eid);
    }
}


#[inline]
fn output_command(com: &Command) {
    match com {
        Command::Stay => {print!("stay\n");},
        Command::Move(v) => {print!("move {}\n", v+1);},
        Command::ChargeFromGrid(c) => {print!("charge_from_grid {}\n", c);},
        Command::ChargeToGrid(c) => {print!("charge_to_grid {}\n", c);},
    }
}

fn main () {
    let timer = Timer::new();

    //graph inputs
    let (v, e) = read!(usize, usize);
    let mut gr = Graph {v: v, e: e, graph: vec![Vec::new(); v], costs: vec![vec![INF/3; v+1]; v+1], next:vec![vec![INF; v]; v]};
    
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
    let fsthalf_minexpect_divt = vec![0; n_pattern];
    let maxexpect_divt = vec![0; n_pattern];
    /*
    let bd = match daytype {
        DayType::Fine | DayType::FineSquall => {10},
        DayType::Rainy | DayType::RainyFine => {20}
    };
    for i in 0..n_pattern {
        let mut mnid = 0;
        let mut mxid = 0;
        for j in 1..n_div+1 {
            if j <= bd && pw_expected[i][mnid] > pw_expected[i][j] {mnid = j;}
            if pw_expected[i][mxid] < pw_expected[i][j] {mxid = j;}
        }
        fsthalf_minexpect_divt[i] = mnid;
        
        if pw_expected[i][mnid] > pw_expected[i][n_div] {
            fsthalf_minexpect_divt[i] = n_div;
        }
        
        maxexpect_divt[i] = mxid;
    }
    */
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
        fsthalf_minexpect_divt: fsthalf_minexpect_divt,
        maxexpect_divt: maxexpect_divt,
    };
    for i in 0..n_grid {
        let (xi, pi) = read!(usize, usize);
        nano_grid.x.push(xi - 1);
        nano_grid.patterns.push(pi - 1);
        nano_grid.pos_to_id[xi - 1] = i;
    }
    
    //EV inputs
    let (n_ev, c_init, c_max, v_max, delta_move) = read!(usize, usize, usize, usize, usize);
    let ev_pos = read!(usize1; n_ev);
    let next: Vec<(usize, usize)> = ev_pos.iter().map(|&x| (x, 0usize)).collect();
    let ev = EVs {
        n_ev: n_ev,
        c_ev_init: c_init, c_ev_max: c_max, v_ev_max: v_max, delta_ev_move: delta_move,
        charge: vec![c_init; n_ev],
        pos: ev_pos,
        next: next,
        adj: vec![Vec::new(); n_ev],
        commands: vec![VecDeque::new(); n_ev],
        pre_commands: vec![Command::Stay; n_ev],
    };

    //other inputs
    let _gamma = read!(f64);
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
        graph: gr,
        nano_grid: nano_grid,
        ev: ev,
        sudden: false,
        complete: false,
        unblanced: false,
        after_sudden: false,
    };

    let bd = match state.nano_grid.day_type {
        DayType::Fine | DayType::FineSquall => {10},
        DayType::Rainy | DayType::RainyFine => {20}
    };
    
    for i in 0..n_pattern {
        let mut mnid = 0;
        let mut mxid = 0;
        for j in 1..n_div+1 {
            if j <= bd && state.nano_grid.pw_expected[i][mnid] > state.nano_grid.pw_expected[i][j] {mnid = j;}
            if state.nano_grid.pw_expected[i][mxid] < state.nano_grid.pw_expected[i][j] {mxid = j;}
        }
        state.nano_grid.fsthalf_minexpect_divt[i] = mnid;
        state.nano_grid.maxexpect_divt[i] = mxid;
    }
    let mut neg_cnt = 0;
    for gid in 0..state.nano_grid.n_grid {
        let pat = state.nano_grid.patterns[gid];
        let mnid = state.nano_grid.fsthalf_minexpect_divt[pat];
        if mnid <= 10 && state.nano_grid.pw_expected[pat][mnid] + (state.nano_grid.c_grid_init as i32) < 0 {neg_cnt += 1;}
    }
    for i in 0..n_pattern {
        if neg_cnt <= 16 || n_ev > 16 {
            if state.nano_grid.pw_expected[i][state.nano_grid.fsthalf_minexpect_divt[i]] > state.nano_grid.pw_expected[i][n_div] {
                state.nano_grid.fsthalf_minexpect_divt[i] = n_div;
            }
        }
    }
    let mut bneed = Vec::new();
    let check_t = match state.nano_grid.day_type {//for one_turn pass
        DayType::Fine | DayType::Rainy => {0usize},
        DayType::FineSquall | DayType::RainyFine => {1usize},
    };
    let mut nxcalc = 0;
    for t in 0..T_MAX {
        read_info(&mut state, t / DIV_RANGE);
        if OUT_INFO && t % DIV_RANGE == check_t {print_status(&state, t);}
        //if t % DIV_RANGE == check_t && ((t / DIV_RANGE) % 4 == 0 || state.sudden) {
        if t % DIV_RANGE == check_t && ((t / DIV_RANGE) == nxcalc || state.sudden) {
            if OUT_INFO {print_status(&state, t);}
            if state.sudden {
                match state.nano_grid.day_type {
                    DayType::RainyFine => {nxcalc = (t / DIV_RANGE) + 1;},
                    _ => {nxcalc = (t / DIV_RANGE) + 4;}
                }
            }
            else {nxcalc = (t / DIV_RANGE) + 4;}
            
            let need_prior = match state.nano_grid.day_type {
                DayType::Fine | DayType::FineSquall => {calc_need_priority_fine(&mut state, (t-check_t) / DIV_RANGE, check_t)},
                DayType::Rainy | DayType::RainyFine => {calc_need_priority_rainy(&mut state, (t-check_t) / DIV_RANGE, check_t)},
            };
            if !state.complete {
                let (perm, mut need) = construct_ev_route(&mut state, need_prior, (t-check_t) / DIV_RANGE);
                construct_commands(&mut state, &perm, &mut need, t, (t-check_t) / DIV_RANGE);
                bneed = need;
            }
        }
        let tt = t;

        for dq in state.nano_grid.chrg_end_time.iter_mut() {
            match dq.front() {
                Some(&ft) => {if tt+1 == ft {dq.pop_front();}},
                None => {}
            }
            match dq.back() {
                Some(&ft) => {if tt+1 == ft {dq.pop_back();}},
                None => {}
            }
        } 
        for (ev, com) in state.ev.commands.iter_mut().enumerate() {
            match com.front() {
                Some(c) => {
                    //let _cc = c.clone();
                    output_command(c);
                    state.ev.pre_commands[ev] = c.clone();
                    com.pop_front();
                    /*
                    match cc {
                        Command::ChargeFromGrid(d) => {
                            match com.front() {
                                Some(Command::ChargeFromGrid(e)) => {},
                                _ => {
                                    if *state.nano_grid.visit_ev_id[state.nano_grid.pos_to_id[d]].front().unwrap() == ev {
                                        state.nano_grid.visit_ev_id[state.nano_grid.pos_to_id[d]].pop_front().unwrap();
                                    }
                                    else {
                                        state.nano_grid.visit_ev_id[state.nano_grid.pos_to_id[d]].pop_back().unwrap();
                                    }
                                }
                            }
                        },
                        Command::ChargeToGrid(d) => {
                            match com.front() {
                                Some(Command::ChargeToGrid(e)) => {},
                                _ => {
                                    if *state.nano_grid.visit_ev_id[state.nano_grid.pos_to_id[d]].front().unwrap() == ev {
                                        state.nano_grid.visit_ev_id[state.nano_grid.pos_to_id[d]].pop_front().unwrap();
                                    }
                                    else {
                                        state.nano_grid.visit_ev_id[state.nano_grid.pos_to_id[d]].pop_back().unwrap();
                                    }
                                }
                            }
                        },
                        _ => {}
                    }
                    */
                },
                None => {
                    output_command(&Command::Stay);
                    state.ev.pre_commands[ev] = Command::Stay.clone();
                }
            }
            
            //print!("stay\n");
        }
        std::io::stdout().flush().unwrap();

        if t >= T_MAX - 1 {break;}
        //if state.complete {continue;}
        if check_t == 1 && t == 0 {continue;}
        
        for i in 0..state.ev.n_ev {
            if state.ev.commands[i].is_empty() {
                construct_additional_commands(&mut state, &mut bneed, i, t, t / DIV_RANGE);
            }
        }
        for i in 0..state.ev.n_ev {
            if state.ev.commands[i].is_empty() {
                construct_additional_commands2(&mut state, &mut bneed, i, t, t / DIV_RANGE);
            }
        }
        /*
        if let DayType::FineSquall = state.nano_grid.day_type {
            for i in 0..state.ev.n_ev {
                if state.ev.commands[i].is_empty() && state.ev.charge[i] > 20000 {
                    let mut mn = 5;
                    let mut nxg = V;
                    for gid in 0..state.nano_grid.n_grid {
                        if state.nano_grid.charge[gid] < 5000 && mn > state.graph.costs[state.ev.next[i].0][state.nano_grid.x[gid]] {
                            mn = state.graph.costs[state.ev.next[i].0][state.nano_grid.x[gid]];
                            nxg = gid;
                        }
                    }
                    if nxg != V && mn != 0 {
                        let mut pos = state.ev.next[i].0; 
                        let acost = match state.ev.pre_commands[i] {
                            Command::Stay => {state.ev.next[i].1},
                            Command::Move(_to) => {
                                if state.ev.next[i].1 == 0 {1} else {state.ev.next[i].1 - 1}
                            },
                            _ => {0}
                        };
                        
                        for _ in 0..acost {
                            state.ev.commands[i].push_back(Command::Move(pos));
                        }
                        let nx = state.nano_grid.x[nxg];
                        while pos != nx {
                            let subnx = state.graph.next[pos][nx];
                            let dist = state.graph.costs[pos][subnx];
                            for _ in 0..dist {
                                state.ev.commands[i].push_back(Command::Move(subnx));
                            }
                            pos = subnx;
                        }
                    }
                }
            }

        }
        */
        
        
    } 
    

    //finally :Tmax info
    read_info(&mut state, 19);

    if OUT_INFO {print_status(&state, 1000);}

    let judge_score = read!(f64);
    if OUT_INFO {eprintln!("{}", judge_score - BASESCORE);}
    
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

