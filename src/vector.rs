#![allow(non_camel_case_types)]

use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone)]
pub struct MVector<f64> {
    pub m: Vec<f64>,
}

impl MVector<f64> {
    fn new(a: Vec<f64>) -> MVector<f64> {
        MVector { m: a }
    }
    // compatability with nalgebra
    pub fn from_row_slice(n: usize, a: &[f64]) -> MVector<f64> {
        let mut m = Vec::new();
        for i in 0..n {
            let c = a[i];
            m.push(c);
        }
        MVector::new(m)
    }

    pub fn as_slice(&self) -> &[f64] {
        self.m.as_slice()
    }

    pub fn len(&self) -> usize {
        self.m.len()
    }
}

// simplified matrix multiplication to avoid using nalgebra
impl Add<MVector<f64>> for MVector<f64> {
    type Output = MVector<f64>;
    fn add(self, rhs: MVector<f64>) -> MVector<f64> {
        let mut ans = Vec::new();
        for (a, b) in self.m.iter().zip(rhs.m) {
            ans.push(a + b);
        }
        MVector::new(ans)
    }
}

impl Add<f64> for MVector<f64> {
    type Output = MVector<f64>;
    fn add(self, rhs: f64) -> MVector<f64> {
        let mut ans = Vec::new();
        for a in self.m {
            ans.push(a + rhs);
        }
        MVector::new(ans)
    }
}

impl Sub<MVector<f64>> for MVector<f64> {
    type Output = MVector<f64>;
    fn sub(self, rhs: MVector<f64>) -> MVector<f64> {
        let mut ans = Vec::new();
        for (a, b) in self.m.iter().zip(rhs.m) {
            ans.push(a - b);
        }
        MVector::new(ans)
    }
}

impl Mul<f64> for MVector<f64> {
    type Output = MVector<f64>;
    fn mul(self, rhs: f64) -> MVector<f64> {
        let mut ans = Vec::new();
        for a in self.m {
            ans.push(a * rhs);
        }
        MVector::new(ans)
    }
}

impl Div<f64> for MVector<f64> {
    type Output = MVector<f64>;
    fn div(self, rhs: f64) -> MVector<f64> {
        let mut ans = Vec::new();
        for a in self.m {
            ans.push(a / rhs);
        }
        MVector::new(ans)
    }
}

impl Mul<MVector<f64>> for f64 {
    type Output = MVector<f64>;
    fn mul(self, rhs: MVector<f64>) -> MVector<f64> {
        let mut ans = Vec::new();
        for a in rhs.m {
            ans.push(a * self);
        }
        MVector::new(ans)
    }
}
