use std::io::Read;

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


#[test]
fn procon_io_test(){
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let mut input = buffer.split_whitespace();

    let n = read!(input, usize);
    for _ in 0..n {
        input!{
            input,
            m: usize,
            mut v: [usize; m],
            c: [chars; m],
        }
        v[m-1] = 99999;
        println!("{}", m);
        println!("{:?}", v);
        println!("{:?}", c);
    }
}