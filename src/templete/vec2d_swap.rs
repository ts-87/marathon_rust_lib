
//mut 1d vec.swap(id1,id2);

pub fn mat_swap<T: Copy>(v: &mut Vec<Vec<T>>, i1: usize, j1: usize, i2: usize, j2: usize) {
    let t = v[i2][j2];
    v[i2][j2] = v[i1][j1];
    v[i1][j1] = t;
}


pub fn mat_swap<T>(v: &mut Vec<Vec<T>>, i1: usize, j1: usize, i2: usize, j2: usize) {
    let p1: *mut T = &mut v[i1][j1];
    let p2: *mut T = &mut v[i2][j2];
    unsafe {
        p1.swap(p2);
    }
}
