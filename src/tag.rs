//! Useful structure for keeping track of the original positions the
//! values come from in a sequence.

use std::cmp::Ordering;
use cycle;

#[derive(Debug, Clone, Copy)]
pub struct Tag {
    pub value: f64,
    pub index: usize,
}

#[allow(dead_code)]
fn remove_cycle_tags(cycles: &[cycle::Cycle<Tag>]) -> Vec<cycle::Cycle<f64>> {
    cycles
        .iter()
        .map(|cycle| cycle::Cycle::new(cycle.max.value, cycle.min.value))
        .collect::<Vec<_>>()
}

impl Tag {
    pub fn new(v: f64, i: usize) -> Tag {
        Tag { value: v, index: i }
    }

    pub fn from(v: &[f64]) -> Vec<Tag> {
        let mut result = Vec::new();
        for (i, value) in v.iter().enumerate() {
            result.push(Tag {
                value: *value,
                index: i,
            });
        }
        result
    }
}

impl PartialEq for Tag {
    fn eq(&self, other: &Tag) -> bool {
        self.value == other.value
    }
}

impl PartialOrd for Tag {
    fn partial_cmp(&self, other: &Tag) -> Option<Ordering> {
        (self.value).partial_cmp(&other.value)
    }
    fn lt(&self, other: &Tag) -> bool {
        self.value < other.value
    }
    fn le(&self, other: &Tag) -> bool {
        self.value <= other.value
    }
    fn gt(&self, other: &Tag) -> bool {
        self.value > other.value
    }
    fn ge(&self, other: &Tag) -> bool {
        self.value >= other.value
    }
}

#[test]
fn test_sort() {
    let t1 = Tag::new(1.0, 0);
    let t2 = Tag::new(-3.0, 1);
    let t3 = Tag::new(5.0, 2);

    let _t = vec![t1, t2, t3];

    assert_eq!(t1 > t2, true)
}
