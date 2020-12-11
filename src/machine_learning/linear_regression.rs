use ml_util::inverse_array2;
use ndarray::{Array1, Array2, ArrayView1, ArrayView2, Axis, NdFloat, arr1, arr2, stack};

use crate::machine_learning::ml_util;

use std::ops::{Add, Sub, Mul};

pub struct LinearRegression<T> 
    where
    T: NdFloat + Add + Sub + Mul
{
    weights: Array1<T>,
}


impl<T: NdFloat + Add + Sub + Mul> LinearRegression<T> {
    pub fn new(n: usize) -> Self {
        Self {weights: Array1::zeros(n)}
    }

    pub fn fit(&mut self, x: &ArrayView2<T>, y: &ArrayView1<T>, intercept: bool) {
        let cx;
        if intercept {
            cx = stack(Axis(1), &[x.view(), Array2::ones((x.shape()[0],1)).view()]).unwrap();  
        }
        else {
            cx = x.add(T::zero());
        }
        let xtx = cx.t().dot(&cx);
        let xtx_inv = inverse_array2(&xtx).unwrap();
        let xty = cx.t().dot(y);
        self.weights = xtx_inv.dot(&xty);
    }

    pub fn fit_grad(&mut self, x: &ArrayView2<T>, y: &ArrayView1<T>, lr: T, epoch: usize, intercept: bool, init: bool) {
        if init {
            let sh = if intercept {x.shape()[1] + 1} else {x.shape()[1]};
            self.weights = Array1::ones(sh);
        }
        let cx;
        if intercept {
            cx = stack(Axis(1), &[x.view(), Array2::ones((x.shape()[0], 1)).view()]).unwrap();
        }
        else {cx = x.add(T::zero());}
        //shaffle...
        for _i in 0..epoch {
            let delta = y - &cx.dot(&self.weights);
            self.weights += &cx.t().dot(&delta).mul(lr);
        }
    }
    
    pub fn fit_ridge(&mut self, x: &ArrayView2<T>, y: &ArrayView1<T>, lambda: T, intercept: bool){
        let cx;
        if intercept {
            cx = stack(Axis(1), &[x.view(), Array2::ones((x.shape()[0],1)).view()]).unwrap();  
        }
        else {
            cx = x.add(T::zero());
        }
        let mut eyep = Array2::eye(cx.shape()[1]);
        eyep[[0,0]] = T::zero();
        let xtx = cx.t().dot(&cx);
        let xtx_inv = inverse_array2(&(&xtx + &eyep.mul(lambda))).unwrap();
        let xty = cx.t().dot(y);
        self.weights = xtx_inv.dot(&xty);
    }

    pub fn predict(&self, x: &ArrayView2<T>, intercept: bool) -> Array1<T> {
        if intercept {
            let cx = stack(Axis(1), &[x.view(), Array2::ones((x.shape()[0], 1)).view()]).unwrap();
            cx.view().dot(&self.weights)
        }
        else {
            x.dot(&self.weights)
        }
    }
}

/*
pub struct LinearRegression_online<T> 
    where
    T: NdFloat,
{
    weights: Array1<T>,
    lr : T
}

impl<T: NdFloat> LinearRegression_online<T> {
    pub fn new(n: usize, lr: T) -> Self {
        Self {weights: Array1::zeros(n), lr: lr}
    }

    pub fn step(&mut self, x: &ArrayView2<T>, y: &ArrayView1<T>, intercept: bool) {

    }
}

impl<T: NdFloat> Regression<T> for LinearRegression_online<T> {
    fn fit(&mut self, x: &ArrayView2<T>, y: &ArrayView1<T>, intercept: bool) {
        if intercept {
            let n = x.shape()[0];
            stack(Axis(1), &[x.view(), Array2::ones((n,1)).view()]).unwrap();
        }
        let xtx = x.t().dot(x);
        let xtx_inv = inverse_array2(&xtx).unwrap();
        let xty =  x.t().dot(y);
        self.weights = xtx_inv.dot(&xty);
    }
    fn predict(&self, x: &ArrayView2<T>, intercept: bool) -> Array1<T> {
        if intercept {
            let n = x.shape()[0];
            stack(Axis(1), &[x.view(), Array2::ones((n,1)).view()]).unwrap();
        }
        x.dot(&self.weights)
    }
}
*/

#[test]
fn test_linear_regression(){
    let x = arr2(&[
        [234.289, 235.6, 159.0, 107.608, 1947., 60.323],
        [259.426, 232.5, 145.6, 108.632, 1948., 61.122],
        [258.054, 368.2, 161.6, 109.773, 1949., 60.171],
        [284.599, 335.1, 165.0, 110.929, 1950., 61.187],
        [328.975, 209.9, 309.9, 112.075, 1951., 63.221],
        [346.999, 193.2, 359.4, 113.270, 1952., 63.639],
        [365.385, 187.0, 354.7, 115.094, 1953., 64.989],
        [363.112, 357.8, 335.0, 116.219, 1954., 63.761],
        [397.469, 290.4, 304.8, 117.388, 1955., 66.019],
        [419.180, 282.2, 285.7, 118.734, 1956., 67.857],
        [442.769, 293.6, 279.8, 120.445, 1957., 68.169],
        [444.546, 468.1, 263.7, 121.950, 1958., 66.513],
        [482.704, 381.3, 255.2, 123.366, 1959., 68.655],
        [502.601, 393.1, 251.4, 125.368, 1960., 69.564],
        [518.173, 480.6, 257.2, 127.852, 1961., 69.331],
        [554.894, 400.7, 282.7, 130.081, 1962., 70.551],
    ]);

    let y = arr1(&[
        83.0, 88.5, 88.2, 89.5, 96.2, 98.1, 99.0, 100.0, 101.2, 104.6, 108.4, 110.8, 112.6,
        114.2, 115.7, 116.9,
    ]);
    
    let mut model = LinearRegression::new(x.shape()[0]);
    model.fit(&x.view(), &y.view(), true);
    let res = model.predict(&x.view(), true);
    println!("{:?}", res);
    


    let mut model = LinearRegression::new(x.shape()[0]);
    model.fit_grad(&x.view(), &y.view(), 1e-8, 20000,true, true);
    let res = model.predict(&x.view(), true);
    println!("{:?}", res);
    

    let mut model = LinearRegression::new(x.shape()[0]);
    model.fit_ridge(&x.view(), &y.view(), 0.02, true);
    let res = model.predict(&x.view(), true);
    println!("{:?}", res);

}
