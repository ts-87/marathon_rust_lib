use ndarray::{Array, Array1, Array2, Axis, NdFloat, arr1, arr2, stack};

use std::ops::{Add, Sub, Mul, Div};
use num_traits::cast::FromPrimitive;

pub struct LogisticRegression<T> 
where
    T: NdFloat + Add + Sub + Mul + Div + FromPrimitive
{
    weights: Array1<T>,
    means: Array1<T>,
    stds: Array1<T>
}


impl<T: NdFloat + Add + Sub + Mul + Div + FromPrimitive> LogisticRegression<T> {
    pub fn new(n: usize) -> Self {
        Self {weights: Array1::zeros(n), means: Array1::zeros(n), stds: Array1::ones(n)}
    }

    pub fn logistic(z: &Array1<T>) -> Array1<T> {
        z.mapv(|x| ((-x).exp() + T::one()).recip())
    }

    pub fn normalize(&mut self, x: &Array2<T>, init: bool) -> Array2<T> {
        if init {
            self.means = x.mean_axis(Axis(0)).unwrap();
            self.stds = x.std_axis(Axis(0), T::one());
        }
        return (x - &self.means).div(&self.stds); 
    }

    fn fit(&mut self, x: &Array2<T>, y: &Array1<T>, lr: T, epoch: usize, normalize: bool, intercept: bool, init: bool) {
        if init {
            let sh = if intercept {x.shape()[1] + 1} else {x.shape()[1]};
            self.weights = Array1::ones(sh);
        }
        let mut cx;
        if normalize {
            cx = self.normalize(x, init);
        }
        else {
            cx = x.add(T::zero());
        }
        if intercept {
            cx = stack(Axis(1), &[cx.view(), Array2::ones((cx.shape()[0],1)).view()]).unwrap();  
            if normalize {
                self.means = stack(Axis(0), &[self.means.view(), Array1::ones(1).view()]).unwrap();
                self.stds = stack(Axis(0), &[self.stds.view(), Array1::ones(1).view()]).unwrap();
            }
        }
        else {
            cx = x.add(T::zero());
        }
        
        

        for _i in 0..epoch {
            let delta = y - &Self::logistic(&cx.dot(&self.weights));
            self.weights += &cx.t().dot(&delta).mul(lr);
        }
    }

    fn predict(&mut self, x: &Array2<T>, normalize: bool, intercept: bool) -> Array1<T> {
        if intercept {
            let mut cx = stack(Axis(1), &[x.view(), Array2::ones((x.shape()[0], 1)).view()]).unwrap();
            if normalize {
                cx = self.normalize(&cx, false);
            }
            Self::logistic(&cx.dot(&self.weights))
        }
        else {
            let mut cx = x.add(T::zero());
            if normalize {
                cx = self.normalize(&cx, false);
            }
            Self::logistic(&cx.dot(&self.weights))
        }
    }
}

use std::fs::File;
use std::io::Read;
use std::path::Path;
#[test]
fn test_logistic_regression(){
    let train_file = Path::new("./src/machine_learning/data/lr_train.data");
    let mut train_data = String::new();

    let mut f = File::open(&train_file).unwrap();
    f.read_to_string(&mut train_data).unwrap();
    let train_iter = train_data.split_whitespace();
    let mut train_x = Vec::new();
    let mut train_y = Vec::new();
    let n = 469;
    let mut cnt = 0;
    for s in train_iter.take(n) {
        let mut it = s.split(',');
        it.next();
        let ch = it.next().unwrap();
        if ch == "B" {
            train_y.push(1.0);
        }
        else {
            train_y.push(0.0);
        }
        it.for_each(|t| train_x.push(t.parse::<f64>().unwrap()));
        cnt += 1;
        if cnt == n {break;}
    }
    let m = train_x.len() / n;
    let mut model = LogisticRegression::new(12);
    model.fit(&Array::from_shape_vec((n,m), train_x).unwrap(), &Array::from_shape_vec(n,train_y).unwrap(), 0.0001, 10000, true, true, true);
    let test_iter = train_data.split_whitespace();
    let mut test_x = Vec::new();
    let mut test_y = Vec::<f64>::new();
    for s in test_iter.skip(n) {
        let mut it = s.split(',');
        it.next();
        let ch = it.next().unwrap();
        if ch == "B" {
            test_y.push(1.0);
        }
        else {
            test_y.push(0.0);
        }
        it.for_each(|t| test_x.push(t.parse::<f64>().unwrap()));
    }
    let n = 100;
    let result = model.predict(&Array::from_shape_vec((n,m), test_x).unwrap(), true, true);
    println!("{:?}", result);
    let mut cnt = 0;
    for (&i, &j) in result.to_vec().iter().zip(test_y.iter()) {
        if i.round() == j.round() {cnt += 1;}
    }
    println!("correct: {}, total: {}, acc: {}", cnt, n, f64::from(cnt) / f64::from(n as i32));
}