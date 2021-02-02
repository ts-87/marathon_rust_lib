
#[inline]
fn modpow (mut x: usize, mut n: usize, m: usize) -> usize {
    let mut ret = 1;
    while n > 0 {
        if (n & 1) == 1 {ret = ret * x % m;}
        x = x * x % m;
        n >>= 1;
    }
    ret
}


pub struct Factrial {
    n: usize,
    m: usize,
    fact: Vec<usize>,
    inv_fact: Vec<usize>,
}
impl Factrial {
    pub fn new(n: usize, m: usize) -> Self {
        let mut fact = vec![0; n+1];
        let mut inv = vec![0; n+1];
        let mut inv_fact = vec![0; n+1];
        inv[1] = 1;
        for i in 2..n+1 {inv[i] = m - (m / i) * inv[m % i] % m;}
        fact[0] = 1;
        inv_fact[0] = 1;
        for i in 1..n+1 {
            fact[i] = fact[i-1] * i % m;
            inv_fact[i] = inv_fact[i-1] * inv[i] % m;
        }
        Self {
            n: n, m: m, fact: fact, inv_fact: inv_fact,
        }
    }

    #[inline]
    pub fn comb (&self, n: usize, k: usize) -> usize {
        if n < k {return 0;}
        self.inv_fact[k] * self.inv_fact[n-k] % self.m * self.fact[n] % self.m
    }
    #[inline]
    pub fn perm (&self, n: usize, k: usize) -> usize {
        if n < k {return 0;}
        self.inv_fact[n-k] * self.fact[n] % self.m
    }
}

#[inline]
fn exgcd(a: i64, b: i64, x: &mut i64, y: &mut i64) -> i64 {
    let mut d = a;
    if b != 0 {
        d = exgcd(b, a%b, y, x);
        *y -= a / b * *x;
    }
    else {
        *x = 1;
        *y = 0;
    }
    d
}
fn mod_inverse(a: i64, m: i64) -> i64 {
    let mut x = 1;
    let mut y = 0;
    exgcd(a, m, &mut x, &mut y);
    (m + x % m) % m
}

#[test]
fn number_test(){

}