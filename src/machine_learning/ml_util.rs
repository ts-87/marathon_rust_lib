use ndarray::{Array1, Array2, arr1, arr2, NdFloat};

pub fn lu_decomposition<T: NdFloat>(a: &Array2<T>, perm: &mut Vec<usize>) -> Option<Array2<T>> {
    let dim = a.shape();
    if dim[0] != dim[1] {return None;}
    let n = dim[0];
    let mut lu = a.clone();
    for i in 0..n-1 {
        let mut mx = lu[[i,i]].abs();
        let mut id = i;
        for j in i+1..n {
            if mx < lu[[j,i]].abs() {
                mx = lu[[j,i]].abs();
                id = j;
            }
        }
        if i != id {
            for j in 0..n {
                lu.swap((i, j),(id, j));
            }
            perm.swap(i, id);
        }
        if !lu[[i,i]].is_normal() {return None;}
        for j in i+1..n {
            let p = lu[[j,i]] / lu[[i,i]];
            for k in i+1..n {
                lu[[j,k]] = lu[[j,k]] - lu[[i,k]] * p;
            }
            lu[[j,i]] = p;
        }
    }
    Some(lu)
}

pub fn solve_linear_equations<T: NdFloat>(a: &Array2<T>, tb: &Array1<T>) -> Option<Array1<T>> {
    if a.shape()[0] != tb.shape()[0] {return None;}

    let n = a.shape()[0];
    let mut perm: Vec<usize> = (0..n).collect();

    let lu;
    match lu_decomposition(a, &mut perm) {
        Some(r) => {lu = r},
        None => {return None;}
    }
    
    let mut b = Array1::zeros(n);
    for i in 0..n {
        b[i] = tb[perm[i]];
    }

    for i in 0..n {
        for j in 0..i {
            b[i] = b[i] - lu[[i,j]] * b[j];
        }
    }

    for i in (0..n).rev() {
        for j in i+1..n {
            b[i] = b[i] - lu[[i,j]] * b[j];
        }
        b[i] /= lu[[i,i]];
    }
    Some(b)
}

pub fn inverse_array2<T: NdFloat>(a: &Array2<T>) -> Option<Array2<T>> {
    let n = a.shape()[0];
    let mut perm: Vec<usize> = (0..n).collect();

    let lu;
    match lu_decomposition(a, &mut perm) {
        Some(r) => {lu = r},
        None => {return None;}
    }
    let n = a.shape()[0];
    let mut b = Array2::zeros((n, n));
    for i in 0..n {
        b[[i, perm[i]]] = T::one();
    }

    for k in 0..n {
        for i in 0..n {
            for j in 0..i {
                b[[i,k]] = b[[i,k]] - lu[[i,j]] * b[[j,k]];
            }
        }

        for i in (0..n).rev() {
            for j in i+1..n {
                b[[i,k]] = b[[i,k]] - lu[[i,j]] * b[[j,k]];
            }
            b[[i,k]] /= lu[[i,i]];
        }
    }
    Some(b)
}


#[test]
fn test_inverse(){
    let a = arr2(
        &[[3., 1., 1., 2.],
        [5., 1., 3., 4.],
        [2., 0., 1., 0.],
        [1., 3., 2., 1.]]);

    match inverse_array2(&a) {
        Some(r) => {
            println!("{:?}", r);
        },
        None => ()
    }

    let a = arr2(
        &[[1., 1., 1., 1.],
        [1., 1., 1., -1.],
        [1., 1., -1., 1.],
        [1., -1., 1., 1.]]);
    let b = arr1(&[0., 4., -4., 2.]);
    match solve_linear_equations(&a, &b) {
        Some(r) => {
            println!("{:?}", r);
        },
        None => ()
    }
    /*
    if solve_linear_equations(&mut a, &b) {
        println!("{:?}", b);
    }
    */
}