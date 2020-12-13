use num_traits::{Float, NumCast, ToPrimitive};
use std::default::Default;

use rand::prelude::*;
use rand::distributions::uniform::SampleUniform;
use rand_pcg::Pcg64Mcg;

#[inline]
pub fn shuffle<T>(v: &mut Vec<T>, rng: &mut Pcg64Mcg) {
    let l = v.len();
    for i in (1..l).rev() {
        let j = rng.gen_range(0, i+1);
        v.swap(i, j);
    }
}
#[derive(Default)]
pub struct Node<T: Float + Default + ToPrimitive +SampleUniform> {
    feature_id: usize,
    split_value: T,
    depth: usize,
    y: usize,
    data: Vec<usize>,
    left: Option<usize>,
    right: Option<usize>
}

#[derive(Default)]
pub struct DecisionTreeClassifier<T> 
where
    T: Float + Default + ToPrimitive + SampleUniform 
{
    nodes: Vec<Node<T>>,
    num_class: usize
}

impl<T: Float + Default + ToPrimitive + SampleUniform> DecisionTreeClassifier<T> {
    pub fn new(num_class: usize) -> Self {
        Self {nodes: Vec::<Node<T>>::default(), num_class: num_class}
    }

    pub fn fit(
        &mut self, x: &Vec<Vec<T>>, y: &Vec<usize>, num_class: usize, max_depth: usize,
        min_features: T, max_features: T, min_samples: T, max_samples: T, seed: u128)
    {
        let mut rng = Pcg64Mcg::new(seed);
        self.num_class = num_class;
        let x_totalsize = x.len();
        let ft_totalsize = x[0].len();
        let size_f_x: T = NumCast::from::<usize>(x_totalsize).unwrap();
        let size_ft_x: T = NumCast::from::<usize>(ft_totalsize).unwrap();
        let r_size = (size_f_x * ((max_samples - min_samples) * rng.gen_range(T::zero(), T::one()) + min_samples)).round();
        let x_size = std::cmp::max(2, NumCast::from::<T>(r_size).unwrap());
        let mut root: Node<T> = Default::default();
        root.depth = 0;
        for _ in 0..x_size {
            root.data.push(rng.gen_range(0, x_totalsize));
        }
        self.nodes.push(root);
        let mut cur_nodeid = 0;

        let mut left_cnt: Vec<usize> = vec![0; num_class];
        let mut right_cnt: Vec<usize> = vec![0; num_class];

        while cur_nodeid < self.nodes.len() {
            let curnode_len = self.nodes.len();
            let curnode = &mut self.nodes[cur_nodeid];
            let mut all_same = true;

            for i in 1..curnode.data.len() {
                if y[curnode.data[i-1]] != y[curnode.data[i]] {
                    all_same = false;
                    break;
                }
            }
            if all_same {
                curnode.y = y[curnode.data[0]];
                cur_nodeid += 1;
                continue;
            }
            if curnode.depth >= max_depth {
                let mut cls_count = vec![0; num_class];
                for &id in curnode.data.iter() {
                    cls_count[id] += 1;
                }
                curnode.y = cls_count.iter().enumerate().max_by(|(_, a), (_, b)| a.cmp(b)).map(|(id, _)| id).unwrap();
                cur_nodeid += 1;
                continue;
            }
            let (mut best_ft_id, mut bestsplitpos) = (num_class, 0);
            let (mut bestscore, mut bestsplitval) = (f64::max_value(), T::zero());
            let r_size = (size_ft_x * ((max_features - min_features) * rng.gen_range(T::zero(), T::one()) + min_features)).round();
            let ft_size = std::cmp::max(2, NumCast::from::<T>(r_size).unwrap());
            let mut feature_ids: Vec<usize> = (0..ft_totalsize).collect();
            shuffle(&mut feature_ids, &mut rng);
            feature_ids.resize(ft_size, 0);
            let mut sorted_ids: Vec<(T, usize)> = Vec::new(); 
            for &fid in feature_ids.iter() {
                sorted_ids.clear();
                for &sid in curnode.data.iter() {
                    sorted_ids.push((x[sid][fid], sid));
                }
                sorted_ids.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let curnode_len = curnode.data.len();
                for i in 0..curnode_len {
                    if i > 0 && sorted_ids[i-1].0 == sorted_ids[i].0 {
                        continue;
                    }
                    for j in 0..num_class {
                        left_cnt[j] = 0;
                        right_cnt[j] = 0;
                    }
                    let (mut cntl, mut cntr) = (0, 0);
                    for j in 0..curnode_len {
                        if j < i {
                            cntl += 1;
                            left_cnt[y[sorted_ids[j].1]] += 1;
                        }
                        else {
                            cntr += 1;
                            right_cnt[y[sorted_ids[j].1]] += 1;
                        }
                    }
                    if cntl == 0 || cntr == 0 {continue;}
                    let mut score = 0.;
                    let mut tmpscore = 1.;
                    for j in 0..num_class {
                        let ratio = left_cnt[j] as f64 / cntl as f64;
                        tmpscore -= ratio * ratio;
                    }
                    score += tmpscore * cntl as f64 / curnode_len as f64;
                    let mut tmpscore = 1.;
                    for j in 0.. num_class {
                        let ratio = right_cnt[j] as f64 / cntr as f64;
                        tmpscore -= ratio * ratio;
                    }
                    score += tmpscore * cntr as f64 / curnode_len as f64;
                    if score < bestscore {
                        bestscore = score;
                        best_ft_id = fid;
                        bestsplitpos = i;
                        if i > 0 {
                            bestsplitval = (sorted_ids[i].0 + sorted_ids[i-1].0) / NumCast::from::<f64>(2.0).unwrap();
                        }
                    }
                }
            }
            if bestsplitpos == 0 {
                let mut cls_count = vec![0; num_class];
                for &id in curnode.data.iter() {
                    cls_count[id] += 1;
                }
                curnode.y = cls_count.iter().enumerate().max_by(|(_, a), (_, b)| a.cmp(b)).map(|(id, _)| id).unwrap();
                cur_nodeid += 1;
                continue;
            }
            
            curnode.feature_id = best_ft_id;
            curnode.split_value = bestsplitval;
            curnode.left = Some(curnode_len);
            curnode.right = Some(curnode_len + 1);
            let mut leftnode: Node<T> = Default::default();
            let mut rightnode: Node<T> = Default::default();
            leftnode.depth = curnode.depth + 1;
            rightnode.depth = curnode.depth + 1;
            for &sid in curnode.data.iter() {
                if x[sid][best_ft_id] < bestsplitval {
                    leftnode.data.push(sid);
                }
                else {
                    rightnode.data.push(sid);
                }
            }
            self.nodes.push(leftnode);
            self.nodes.push(rightnode);
            cur_nodeid += 1;
        }
    }

    pub fn predict(&self, x: &Vec<T>) -> usize {
        let mut curnode_id = Some(0);
        let mut target_id = 0;
        while let Some(nxnode_id) = curnode_id {
            target_id = nxnode_id;
            if x[self.nodes[nxnode_id].feature_id] < self.nodes[nxnode_id].split_value {
                curnode_id = self.nodes[nxnode_id].left;
            }
            else {
                curnode_id = self.nodes[nxnode_id].right;
            }
        }
        self.nodes[target_id].y
    }
}

#[derive(Default)]
pub struct RandomForestClassifier<T> 
where
    T: Float + Default + ToPrimitive + SampleUniform
{
    trees: Vec<DecisionTreeClassifier<T>>,
    num_class: usize,
    num_trees: usize
}

impl<T: Float + Default + ToPrimitive + SampleUniform> RandomForestClassifier<T> {
    pub fn new(num_class: usize) -> Self {
        Self {trees: Vec::<DecisionTreeClassifier<T>>::new(), num_class: num_class, num_trees: 0}
    }

    pub fn fit(&mut self, x: &Vec<Vec<T>>, y: &Vec<usize>, num_trees: usize, num_class: usize,
        max_depth: usize, min_features: T, max_features: T, min_samples: T, max_samples: T) {
        self.num_class = num_class;
        self.num_trees = num_trees;
        for i in 0..num_trees {
            let mut tree = DecisionTreeClassifier::<T>::new(num_class);
            tree.fit(x, y, self.num_class, max_depth, min_features, max_features, min_samples, max_samples, i as u128);
            self.trees.push(tree);
        }
    }
    
    pub fn predict(&self, x: &Vec<T>) -> usize {
        let mut cls_count = vec![0; self.num_class];
        for i in 0..self.num_trees {
            cls_count[self.trees[i].predict(x)] += 1;
        }
        cls_count.iter().enumerate().max_by(|(_, a), (_, b)| a.cmp(b)).map(|(id, _)| id).unwrap()
    }
    /*
    pub fn save(&self) {
        println!("{}" , self.trees.len());
        for i in 0..self.num_trees {
            println!("{}", self.trees[i].nodes.len());
            for node in self.trees[i].nodes.iter() {
                println!
                ("{} {} {:?} {:?} {:?} {}",
                node.feature_id,
                node.depth,
                node.split_value,
                node.left,
                node.right,
                node.y)
            }
        }
    }
    */
}


use std::fs::File;
use std::io::Read;
use std::path::Path;
#[test]
fn test_decision_tree(){
    let train_file = Path::new("./src/machine_learning/data/dtcls_train");
    let mut train_data = String::new();
    let mut f = File::open(&train_file).unwrap();
    f.read_to_string(&mut train_data).unwrap();
    let train_iter = train_data.split_whitespace();

    let mut tmp_x = Vec::new();
    let mut tmp_y = Vec::new();
    let n = 150;
    let m = 4;
    for s in train_iter {
        let mut it = s.split(',');
        let mut xi = Vec::new(); 
        for _ in 0..m {
            let val: f64 = it.next().unwrap().parse().unwrap();
            xi.push(val);
        }
        tmp_x.push(xi);
        let yi = match it.next().unwrap() {
            "Iris-setosa" => 0usize,
            "Iris-versicolor" => 1,
            _ => 2
        };
        tmp_y.push(yi);
    }

    let mut perm: Vec<usize> = (0..n).collect();
    shuffle(&mut perm, &mut Pcg64Mcg::new(42));
    let (num_train, num_test) = (100, n - 100);
    let mut train_x = Vec::new();
    let mut train_y = Vec::new();
    let mut test_x = Vec::new();
    let mut test_y = Vec::new();
    for i in 0..num_train {
        train_x.push(tmp_x[perm[i]].clone());
        train_y.push(tmp_y[perm[i]]);
    }
    for i in num_train..num_test + num_train {
        test_x.push(tmp_x[perm[i]].clone());
        test_y.push(tmp_y[perm[i]]);
    }
    /*
    let mut model = DecisionTreeClassifier::new(m);
    model.fit(&train_x, &train_y, 4, 100, 0.5, 1.0, 1.0, 1.0, 42);
    */
    let mut model = RandomForestClassifier::<f64>::new(m);
    model.fit(&train_x, &train_y, 30, 4, 10, 0.5, 1.0, 0.7, 1.0);
    let mut correct_cnt = 0;
    for i in 0..num_test {
        if test_y[i] == model.predict(&test_x[i]) {
            correct_cnt += 1;
        }
    }
    println!("{}", correct_cnt as f64 / num_test as f64);
}
