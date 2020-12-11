pub struct DisjointSet {
    numtree: usize,
    data: Vec<usize>,
    size: Vec<usize>
}

impl DisjointSet {
    pub fn new(n: usize) -> Self{
        Self {
            numtree: n,
            data: (0..n).map(|i| i).collect::<Vec<usize>>(),
            size: vec![1; n],
        }
    }

    pub fn find(&mut self, x: usize) -> usize {
        if self.data[x] == x {
            x
        }
        else {
            let px = self.data[x];
            self.data[x] = self.find(px);
            self.data[x]
        }
    }

    pub fn unite(&mut self, x: usize, y: usize) {
        let mut px = self.find(x);
        let mut py = self.find(y);
        if px == py {
            return;
        }
        if self.size[py] < self.size[px] {
            std::mem::swap(&mut px, &mut py);
        }
        self.data[px] = py;
        self.size[py] += self.size[px];
        self.numtree -= 1;
    }

    pub fn same(&mut self, x: usize, y: usize) -> bool {
        let px = self.find(x);
        let py = self.find(y);
        self.data[px] == self.data[py]
    }

    pub fn tree_size(&mut self, x: usize) -> usize {
        let px = self.find(self.data[x]);
        self.size[px]
    }

}

#[test]
fn test_disjoint_set(){
    let mut uf = DisjointSet::new(11);
    let inpt: [(usize,usize); 10] = [(0,1),(4,5),(6,8),(9,10),(1,7),(3,6),(7,9),(1,10),(3,4),(4,5)];
    for (i,j) in inpt.iter() {
        uf.unite(*i,*j);
    }
    for i in 0..11 {
        for j in 0..11 {
            println!("{}-{} :{}",i,j,uf.same(i,j));
        }
    }
    for i in 0..11 {
        println!("size{}: {}",i,uf.tree_size(i));
    }
    println!("num_tree: {}",uf.numtree);
}

use std::io::Read;
#[test]
fn test_kruskal() {
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let mut iter = buffer.split_whitespace();

    let n: usize = iter.next().unwrap().parse().unwrap();
    let e: usize = iter.next().unwrap().parse().unwrap();

    let mut edge: Vec<(usize, usize, i32)> = (0..e).map(|_| {
        (iter.next().unwrap().parse().unwrap(),
        iter.next().unwrap().parse().unwrap(),
        iter.next().unwrap().parse().unwrap())
    }).collect();

    let mut uf = DisjointSet::new(n);
    let mut tot = 0;
    edge.sort_by_key(|e| e.2);
    for (sr, ds, w) in edge.iter() {
        if !uf.same(*sr, *ds) {
            uf.unite(*sr, *ds);
            tot += w;
        }
    }
    println!("{}", tot);
}