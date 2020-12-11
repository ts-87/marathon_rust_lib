pub trait BinarySearch<T> {
    fn lower_bound(&self, x: &T) -> usize;
    fn upper_bound(&self, x: &T) -> usize;
}

impl<T: Ord> BinarySearch<T> for [T] {
    fn lower_bound(&self, x: &T) -> usize {
        let mut lo = 0;
        let mut hi = self.len();

        while lo != hi {
            let mid = (lo + hi) >> 1;
            if self[mid] < *x {
                lo = mid + 1;
            }
            else {
                hi = mid;
            }
        }
        lo
    }

    fn upper_bound(&self, x: &T) -> usize {
        let mut lo = 0;
        let mut hi = self.len();

        while lo != hi {
            let mid = (lo + hi) >> 1;
            if self[mid] > *x {
                hi = mid;
            }
            else {
                lo = mid + 1;
            }
        }
        lo
    }
}

#[test]
fn test_binarysearch(){
    let vec = vec![1, 2, 4, 6, 7, 12, 54, 60];
    assert_eq!(vec.lower_bound(&4), 2);
    assert_eq!(vec.upper_bound(&4), 3);

    let arr = [3, 7, 8, 11, 15, 22, 26];
    assert_eq!(arr.lower_bound(&11), 3);
    assert_eq!(arr.upper_bound(&11), 4);
}