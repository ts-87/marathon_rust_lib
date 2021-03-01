#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_macros)]

use std::{io::Read, vec};
macro_rules! input {
    ($iter:expr) => {};
    ($iter:expr, ) => {};

    ($iter:expr, $var:ident : $t:tt $($r:tt)*) => {
        let $var = read!($iter, $t);
        input!{$iter $($r)*}
    };

    ($iter:expr, mut $var:ident : $t:tt $($r:tt)*) => {
        let mut $var = read!($iter, $t);
        input!{$iter $($r)*}
    };
}

macro_rules! read {
    ($iter:expr, ( $($t:tt),* )) => {
        ( $(read!($iter, $t)),* )
    };

    ($iter:expr, [ $t:tt ; $len:expr ]) => {
        (0..$len).map(|_| read!($iter, $t)).collect::<Vec<_>>()
    };

    ($iter:expr, chars) => {
        read!($iter, String).chars().collect::<Vec<char>>()
    };

    ($iter:expr, usize1) => {
        read!($iter, usize) - 1
    };

    ($iter:expr, $t:ty) => {
        $iter.next().unwrap().parse::<$t>().expect("Parse error")
    };
}

use std::cmp::Ordering;
struct State {
    dist: usize,
    cost: usize,
    p0: usize,
    g: Vec<usize>,
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some((other.cost + other.dist).cmp(&(self.cost + self.dist)))
    }
}
impl Eq for State {}
impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        (other.cost + other.dist).eq(&(self.cost + self.dist))
    }
}
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.cost + other.dist).cmp(&(self.cost + self.dist))
    }
}
use std::collections::{BinaryHeap, HashMap};
const INF: usize = 100000;

//https://atcoder.jp/contests/indeednow-quala/tasks/indeednow_2015_quala_4
fn main() {
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let mut input = buffer.split_whitespace();
    
    input!{
        input,
        h: usize,
        w: usize,
    }
    let mut g = vec![0; h*w];
    let mut p0 = 0;
    for i in 0..h*w {
        let ci = read!(input, usize);
        g[i] = ci;
        if ci == 0 {p0 = i;}
    }
    let mut dist = HashMap::new();
    let mut pq = BinaryHeap::new();
    dist.insert(g.clone(), 0);
    let cost = g.iter().enumerate().filter(|(_, c)| **c != 0)
            .map(|(i, &c)| 
            ((i/w)as i32 - ((c-1)/w) as i32).abs() + ((i%w) as i32 - ((c-1)%w) as i32).abs())
            .sum::<i32>() as usize;
    //eprintln!("{}", cost);
    if cost == 0 {
        print!("0\n");
        return;
    }
    pq.push(State {dist: 0, cost: cost, p0: p0, g: g});
    let mut ans = 24;
    while let Some(mut st) = pq.pop() {
        let r = st.p0 / w;
        let c = st.p0 % w;
        for d in [!0, 0, 1, 0, !0].windows(2) {
            let nr = r + d[0];
            let nc = c + d[1];
            if nr >= h || nc >= w {continue;}
            st.g.swap(r*w + c, nr*w + nc);
            let co = if let Some(co) = dist.get(&st.g) {*co} else {INF};
            if st.dist + 1 < co {
                let num = (st.g[r*w+c] - 1) as i32;
                let diff = (r as i32 - num / w as i32).abs() + (c as i32 - num % w as i32).abs()
                        - (nr as i32 - num / w as i32).abs() - (nc as i32 - num % w as i32).abs();
                dist.insert(st.g.clone(), st.dist + 1);
                let totcost = st.cost as i32 + diff;
                if totcost == 0 {
                    ans = ans.min(st.dist + 1);
                    print!("{}\n", ans);
                    return;
                }
                pq.push(State {dist: st.dist + 1, cost: totcost as usize, p0: nr*w+nc, g: st.g.clone()});
            }
            st.g.swap(r*w + c, nr*w + nc);
        }
    }
}

