// NonNan numbers Floating point numbers cannot be sorted because they
// can be Nan.  We introduce a new type of number guaranteed not to be
// an NaN to allow floats to be sorted or the maximum found.

use std::cmp::Ordering;

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct NonNan(f64);

impl NonNan {
    pub fn new(val: f64) -> Option<NonNan> {
        if val.is_nan() {
            None
        } else {
            Some(NonNan(val))
        }
    }
}

impl Eq for NonNan {}
impl Ord for NonNan {
    fn cmp(&self, other: &NonNan) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
