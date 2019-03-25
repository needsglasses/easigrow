/// calculates the combination type product of a range of variables
/// [a b c], 3 gives aaa, baa caa aba bba cba aca bca cca aac ....
/// where np is the size of the combination tupple
pub fn permutations<T: Copy>(range: &[T], np: usize) -> Vec<Vec<T>> {
    let n = range.len();
    let mut factors = Vec::new();

    for i in 0..n.pow(np as u32) {
        let mut base = Vec::new();
        for j in 0..np as u32 {
            base.push(range[((i / n.pow(j)) % n) as usize]);
        }
        factors.push(base);
    }

    factors
}

#[test]
fn test_permutation() {
    let x = vec![1.0f64, 2.0, 3.0];
    let ans = vec![
        vec![1.0f64, 1.0],
        vec![2.0, 1.0],
        vec![3.0, 1.0],
        vec![1.0, 2.0],
        vec![2.0, 2.0],
        vec![3.0, 2.0],
        vec![1.0, 3.0],
        vec![2.0, 3.0],
        vec![3.0, 3.0],
    ];

    assert!(permutations(&x, 2).eq(&ans));
}
