use ndarray::{Array, Array1, Array2, Array3, Axis, NdFloat, Dimension};
use std::ops::{Add, Sub, Mul, Div};
use num_traits::{FromPrimitive};


pub struct Normalize<T> 
where
    T: NdFloat + Add + Sub + Mul + Div + FromPrimitive
{
    means: Array1<T>,
    stds: Array1<T>
}
impl<T: NdFloat + Add + Sub + Mul + Div + FromPrimitive> Normalize<T> {
    pub fn new(m: usize) -> Self{
        Self {
            means: Array1::<T>::zeros(m),
            stds: Array1::<T>::ones(m)
        }
    }

    pub fn normalize(&mut self, x: &Array2<T>, init: bool) -> Array2<T> {
        if init {
            self.means = x.mean_axis(Axis(0)).unwrap();
            self.stds = x.std_axis(Axis(0), T::one());
        }
        return (x - &self.means).div(&self.stds); 
    }
}

pub fn relu<T: NdFloat, D: Dimension>(x: &Array<T, D>) -> Array<T, D> {
    x.mapv(|a| if a > T::zero() {a} else {T::zero()})
}

pub fn sigmoid<T: NdFloat, D: Dimension>(x: &Array<T, D>) -> Array<T, D> {
    x.mapv(|a| ((-a).exp() + T::one()).recip())
}

//regression
pub struct NeuralNetwork<T> 
where
    T: NdFloat + Add + Sub + Mul + Div
{
    w1: Array2<T>,
    c1: Array2<T>,
    w2: Array2<T>,
    c2: Array2<T>
}

impl<T: NdFloat + Add + Sub + Mul + Div> NeuralNetwork<T> {

    pub fn new() -> Self{
        Self {
            w1: Array2::<T>::zeros((1, 1)),
            c1: Array2::<T>::zeros((1, 1)),
            w2: Array2::<T>::zeros((1, 1)),
            c2: Array2::<T>::zeros((1, 1))
        }
    }

    pub fn fit(&mut self, x: &Array2<T>, y: &Array2<T>, hidden_dim: usize, lr: T, epoch: usize, _seed: u128) {
        let n = x.shape()[0];
        let m = x.shape()[1];
        let ydim = y.shape()[1];

        let init_f: T = T::from(0.03).unwrap();
        self.w1 = Array2::<T>::from_elem((hidden_dim, m), init_f);
        self.c1 = Array2::<T>::from_elem((hidden_dim, 1),init_f);
        self.w2 = Array2::<T>::from_elem((ydim, hidden_dim), init_f);
        self.c2 = Array2::<T>::from_elem((ydim, 1), init_f);

        let mut h1 = self.w1.dot(&x.t()).add(&self.c1);
        let mut z1 = relu(&h1);
        let mut h2 = self.w2.dot(&z1).add(&self.c2);
        let mut yhat = sigmoid(&h2);

        for _e in 0..epoch {
            let mut dl_dw2 = Array2::<T>::zeros((ydim, hidden_dim));
            let mut dl_dc2 = Array2::<T>::zeros((1, ydim));
            let mut dl_dw1 = Array2::<T>::zeros((hidden_dim, m));
            let mut dl_dc1 = Array2::<T>::zeros((1, hidden_dim));
            //let f2: T = T::from(2.0).unwrap();
            let neg1: T = T::from(-1.0).unwrap();
            for i in 0..n {
                //let dl_dy = yhat.slice_mut(s![.., i..i+1]).t().sub(&y.slice(s![i..i+1, ..])).mul(f2);
                let mut dl_dy = -(y.slice(s![i..i+1, ..])).div(&yhat.slice(s![.., i..i+1]).t());
                dl_dy += &(y.slice(s![i..i+1, ..]).mul(neg1).add(T::one())).div(&(yhat.slice(s![.., i..i+1]).mul(neg1).add(T::one())).t());

                //let dy_dh2 = Array2::<T>::eye(ydim);
                let dy_dh2 = Array2::<T>::eye(ydim).mul(&sigmoid(&h2.slice(s![.., i..i+1]).add(T::zero())).mul(sigmoid(&h2.slice(s![.., i..i+1]).mul(neg1).add(T::one()))));
                let dh2_dc2 = Array2::<T>::eye(ydim);
                let mut dh2_dw2 = Array3::<T>::zeros((ydim, ydim, hidden_dim));
                for j in 0..ydim {
                    dh2_dw2.slice_mut(s![j, j, ..]).assign(&z1.slice(s![.., i]).t());
                }
                let dh2_dz2 = self.w2.clone();

                let dz1_dh1 = Array2::<T>::eye(hidden_dim).mul(&h1.slice(s![.., i..i+1]).mapv(|x| if x > T::zero() {T::one()} else {T::zero()}));
                let dh1_dc1 = Array2::<T>::eye(hidden_dim);
                let mut dh1_dw1 = Array3::<T>::zeros((hidden_dim, hidden_dim, m));
                for j in 0..hidden_dim {
                    dh1_dw1.slice_mut(s![j, j, ..]).assign(&x.slice(s![i, ..]));
                }

                let dl_dh2 = dl_dy.dot(&dy_dh2);
                for k in 0..ydim {
                    let tmp = &dl_dw2.slice(s![k..k+1, ..]).add(&dl_dh2.dot(&dh2_dw2.slice(s![k, .., ..])));
                    dl_dw2.slice_mut(s![k..k+1, ..]).assign(tmp);
                }
                //dl_dw2 = dl_dh2.dot(&dh2_dw2);
                dl_dc2 += &dl_dh2.dot(&dh2_dc2);
                let dl_dh1 = dl_dh2.dot(&dh2_dz2).dot(&dz1_dh1);
                for k in 0..hidden_dim {
                    let tmp = &dl_dw1.slice(s![k..k+1, ..]).add(&dl_dh1.dot(&dh1_dw1.slice(s![k, .., ..])));
                    dl_dw1.slice_mut(s![k..k+1, ..]).assign(tmp);
                }
                //dl_dw1 = dl_dh1.dot(&dh1_dw1);
                dl_dc1 += &dl_dh1.dot(&dh1_dc1);
            }
            self.w1 -= &dl_dw1.mul(lr);
            self.c1 -= &dl_dc1.into_shape((hidden_dim, 1)).unwrap().mul(lr);
            self.w2 -= &dl_dw2.mul(lr);
            self.c2 -= &dl_dc2.into_shape((ydim, 1)).unwrap().mul(lr);

            h1 = self.w1.dot(&x.t()).add(&self.c1);
            z1 = relu(&h1);
            h2 = self.w2.dot(&z1).add(&self.c2);
            yhat = sigmoid(&h2);
        }
    }

    pub fn predict(&self, x: &Array2<T>) -> Array2<T> {
        let h1 = &self.w1.dot(&x.t()).add(&self.c1);
        let z1 = &relu(h1);
        let h2 = &self.w2.dot(z1).add(&self.c2);
        sigmoid(h2).reversed_axes()
    }
}

use::ndarray::{arr1, arr2};
use std::fs::File;
use std::io::Read;
use std::path::Path;
#[test]
fn test_neural_network(){
    let x = arr1(&[1.,2.,-1.,4.,-6.]);
    println!("{:?}", relu(&x));
    println!("{:?}", sigmoid(&x));
    let x = arr2(&[[1.,2.,-1.,4.,-6.],[1.,-2.,3.,-4.,5.]]);
    println!("{:?}", relu(&x));
    println!("{:?}", sigmoid(&x));

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
            train_y.push(0.0);
        }
        else {
            train_y.push(0.0);
            train_y.push(1.0);
        }
        it.for_each(|t| train_x.push(t.parse::<f64>().unwrap()));
        cnt += 1;
        if cnt == n {break;}
    }
    let m = train_x.len() / n;

    let mut normalizer = Normalize::<f64>::new(m);
    let train_x_arr = normalizer.normalize(&Array::from_shape_vec((n,m), train_x).unwrap(), true);
    
    let mut model = NeuralNetwork::new();
    model.fit(&train_x_arr, &Array::from_shape_vec((n, 2),train_y).unwrap(), 10, 0.00001, 1000, 42);
    
    let test_iter = train_data.split_whitespace();
    let mut test_x = Vec::new();
    let mut test_y = Vec::<(f64, f64)>::new();
    for s in test_iter.skip(n) {
        let mut it = s.split(',');
        it.next();
        let ch = it.next().unwrap();
        if ch == "B" {
            test_y.push((1.0, 0.0));
        }
        else {
            test_y.push((0.0, 1.0));
        }
        it.for_each(|t| test_x.push(t.parse::<f64>().unwrap()));
    }
    let n = 100;

    let test_x_arr = normalizer.normalize(&Array::from_shape_vec((n,m), test_x).unwrap(), false);
    
    let result = model.predict(&test_x_arr);
    println!("{:?}", result);
    let mut cnt = 0;
    for i in 0..n {
        if result[[i,0]] < result[[i,1]] && test_y[i].0 < test_y[i].1 || result[[i,0]] > result[[i,1]] && test_y[i].0 > test_y[i].1 {
            cnt += 1;
        }
    }
    println!("correct: {}, total: {}, acc: {}", cnt, n, f64::from(cnt) / f64::from(n as i32));
}