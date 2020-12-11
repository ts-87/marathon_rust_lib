use std::ops::{Add, AddAssign};

pub struct FenwickTree<T>
where
    T: AddAssign + Add<Output = T> + Copy,
{
    n: usize,
    data: Vec<T>,
    e: T
}

impl<T> FenwickTree<T>
where
    T: AddAssign + Add<Output = T> + Copy, 
{
    pub fn new(n: usize, e: T) -> Self{
        Self {
            n: n,
            data: vec![e; n+1],
            e: e,
        }
    }

    pub fn add(&mut self, mut id: usize, x: T) {
        while id <= self.n {
            self.data[id] += x;
            id += id & !id + 1;
        }
    }

    pub fn sum(&self, mut id: usize) -> T {
        let mut res: T = self.e;
        while id != 0 {
            res += self.data[id];
            id &= id - 1;
        }
        res
    }
}

use std::io::Read;
#[test]
fn test_fenwick_tree(){
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let mut iter = buffer.split_whitespace();

    let n: usize = iter.next().unwrap().parse().unwrap();
    let q: usize = iter.next().unwrap().parse().unwrap();
    let mut ft = FenwickTree::new(n, 0);

    for _ in 0..q {
        let (com, x, y): (usize, usize, usize) = { 
            (iter.next().unwrap().parse().unwrap(),
            iter.next().unwrap().parse().unwrap(),
            iter.next().unwrap().parse().unwrap())
        };
        if com == 0 {
            ft.add(x, y);
        }
        else {
            let ret = ft.sum(y) - ft.sum(x-1);
            print!("{}\n", ret);
        }
    }
}