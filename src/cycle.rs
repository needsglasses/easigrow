//! Collection of fatigue cycle counting routines such as rainflow,
//! turning point and filtering.

/// Paul White (Dec 2014)
use tag;
use std::ops::Sub;
use std;
use std::iter::Iterator;

use log::info;

#[derive(Debug, Clone)]
pub enum CycleMethod {
    Rainflow,
    Tension,
}

#[derive(Debug, Clone)]
pub struct Limit {
    pub max: f64,
    pub min: f64,
}

impl Limit {
    pub fn from_vec(limits: &[f64]) -> Limit {
        Limit {
            max: limits[0].max(limits[1]),
            min: limits[0].min(limits[1]),
        }
    }

    pub fn new(a: f64, b: f64) -> Limit {
        if a > b {
            Limit { max: a, min: b }
        } else {
            Limit { max: b, min: a }
        }
    }
}

/// Data for the changes or limitations applied to the sequence.
#[derive(Debug, Clone)]
pub struct SequenceModifiers {
    /// limit any of the turning points to a maximum of the cap value.
    pub cap_max: Option<f64>,
    /// limit any of the turning points to a minimum of the cap value.
    pub cap_min: Option<f64>,
    /// remove any of the turning points larger than this value.
    pub remove_bigger: Option<f64>,
    /// remove any of the turning points smaller than this value.
    pub remove_smaller: Option<f64>,
    /// If true only keep turning points for remaining cycles
    pub cycles: bool,
    /// Close-up all of the rainflow cycles.
    pub reorder: bool,
    /// Remove non-turning points.
    pub turning_points: bool,
    /// Name of file to write modified values too.
    pub outfile: Option<String>,
}

/// Data for the changes or limitations applied to the cycles.
#[derive(Debug, Clone)]
pub struct CycleModifiers {
    /// cap the cycles to a maximum range of this value.
    pub cap_max: Option<f64>,
    /// cap the cycles to a minimum range of this value.
    pub cap_min: Option<f64>,
    /// Apply a rise fall criteria, cycles smaller than this will be removed.
    pub remove_smaller: Option<f64>,
    /// Remove cycles that are bigger than this.
    pub remove_bigger: Option<f64>,
    /// Remove cycles within this region.
    pub remove_region: Option<Limit>,
    /// Name of file to write modified values too.
    pub outfile: Option<String>,
}

#[derive(Debug, Copy, Clone)]
pub struct Cycle<T> {
    pub max: T,
    pub min: T,
}

/// filter out smaller cycles
impl<T: PartialOrd + Sub<Output = T>> Cycle<T> {
    pub fn range(self) -> T {
        self.max - self.min
    }
    pub fn within(&self, upper: T, lower: T) -> bool {
        (self.max < upper) & (self.min > lower)
    }
}

impl<T: PartialOrd> Cycle<T> {
    pub fn new(a: T, b: T) -> Cycle<T> {
        if a > b {
            Cycle { max: a, min: b }
        } else {
            Cycle { max: b, min: a }
        }
    }
}

impl Cycle<tag::Tag> {
    /// Finds the range of a cycle.
    pub fn range(self) -> f64 {
        self.max.value - self.min.value
    }
    pub fn r(self) -> f64 {
        self.min.value / self.max.value
    }
    /// Tests if the cycle lies between upper and lower limits.
    pub fn within(&self, upper: f64, lower: f64) -> bool {
        (self.max.value <= upper) & (self.min.value >= lower)
    }
}

// Collection of filters to apply to the sequence prior to cycle counting.
pub fn process_seq_mods(sequence: &[tag::Tag], seq_mods: &SequenceModifiers) -> Vec<tag::Tag> {
    let mut mod_seq = if seq_mods.reorder {
        reorder_sequence(sequence)
    } else {
        sequence.to_vec()
    };

    if seq_mods.turning_points {
        mod_seq = turning_points(&mod_seq);
    };

    mod_seq.iter()
    // cap maximum value
        .map(|&s| if let Some(max) = seq_mods.cap_max {
            if s.value > max {
                tag::Tag { value: max, index: s.index }
            } else { s }
        } else { s })
        
    // cap minimum value
        .map(|s| if let Some(min) = seq_mods.cap_min {
            if s.value < min {
                tag::Tag { value: min, index: s.index }
            } else { s }
        } else { s })
        
    // remove smaller points
        .filter(|s| if let Some(min) = seq_mods.remove_smaller {
            s.value < min 
        } else { true })
        
    // remove bigger points
        .filter(|s| if let Some(max) = seq_mods.remove_bigger {
            s.value > max
        } else { true })
        
        .collect::<Vec<_>>()
}

// Collection of filters to apply to the list of cycles prior to using for crack growth.
pub fn process_cycle_mods(
    cycles: &[Cycle<tag::Tag>],
    cycle_mods: &CycleModifiers,
) -> Vec<Cycle<tag::Tag>> {
    cycles.iter()
    // risefall
        .filter(|c| if let Some(ref risefall) = cycle_mods.remove_smaller {
            c.range() > *risefall
        } else { true })

        .filter(|c| if let Some(ref bigger) = cycle_mods.remove_bigger {
            c.range() < *bigger
        } else { true })

    // deadband
        .filter(|c| if let Some(ref deadband) = cycle_mods.remove_region {
            ! (c.max.value <= deadband.max && c.min.value >= deadband.min)
        } else { true })

        .map(|&c| if let Some(max) = cycle_mods.cap_max {
            if c.max.value > max {
                Cycle { max: tag::Tag { value: max, index: c.max.index }, min: c.min }
            } else { c }
        } else { c })

        .map(|c| if let Some(min) = cycle_mods.cap_min {
            if c.min.value < min {
                Cycle { max: c.max, min: tag::Tag { value: min, index: c.min.index } }
            } else { c }
        } else { c })

        .collect::<Vec<_>>()
}

/// Re-create a sequence from a vector of cycles.
pub fn sequence_from_cycles(cycles: &[Cycle<tag::Tag>]) -> Vec<tag::Tag> {
    // find the maximum possible index value
    let max_index = cycles.iter().fold(0.0, |m, c| {
        (c.max.index as f32).max(c.min.index as f32).max(m)
    }) as usize;
    let mut seq = Vec::with_capacity(max_index + 1);

    // initialise the sequence
    for _ in 0..=max_index {
        seq.push(tag::Tag::new(0.0, max_index))
    }

    // assign each cycle to the turning points
    for cycle in cycles {
        seq[cycle.min.index].index = cycle.min.index;
        seq[cycle.min.index].value = cycle.min.value;

        seq[cycle.max.index].index = cycle.max.index;
        seq[cycle.max.index].value = cycle.max.value;
    }

    // filter out any turning points with unassigned points
    seq.into_iter()
        .filter(|s| s.index >= max_index)
        .collect::<Vec<tag::Tag>>()
}

/// Extracts out the turning points from a sequence.
///
/// The sequence may be any type that can be sorted. As such we can
/// use a tagged structure to keep the original position of the
/// turning point in the sequence vector.
pub fn turning_points<T: PartialOrd + Clone>(seq: &[T]) -> Vec<T> {
    let mut tp = Vec::new();
    let mut s = seq.iter();

    let mut a = match s.next() {
        Some(x) => x,
        None => return vec![],
    };
    let mut b = match s.next() {
        Some(x) => x,
        None => return vec![a.clone()],
    };
    tp.push(a.clone());

    loop {
        let c = match s.next() {
            Some(x) => x,
            None => {
                tp.push(b.clone());
                break;
            }
        };

        if !between_eq(a, b, c) {
            tp.push(b.clone());
            a = b;
        }
        b = c;
    }
    tp
}

#[inline]
fn between_eq<T: PartialOrd>(a: &T, b: &T, c: &T) -> bool {
    (a <= b && b <= c) || (a >= b && b >= c)
}

/// Rainflow cycle counting.
///
/// This routine takes a turning point sequence `x` and returns a
/// tuple of
/// (the max and min of each closed cycle, remaining unclosed turning points).
///
/// A closed rainflow cycle is any adjacent pair of turning points b,c that
/// lie between turning points a and d
///```ignore
///            d
///         b  /
///      \  /\/
///       \/  c
///        a
///```
/// This routine repeatedly extracts pairs of adjacent cycles that
/// are an interruption or rainflow cycle from the sequence until no
/// more can be extracted. This may leave a sequence of leftover
/// points which are increasing and/or decreasing set of cycles just
/// like a wavelet since no pairs of points lie in between adjacent
/// points. In order to close up all cycles the sequence should
/// firstly be reordered.
///
/// There seems to be no difference between range pair counting and
/// rainflow counting. The same definition of a cycle as an
/// intermediate cycle that lies between two turning points is a range
/// pair or rainflow cycle.
pub fn rainflow(seq: &[tag::Tag]) -> (Vec<Cycle<tag::Tag>>, Vec<tag::Tag>) {
    let (mut cycles, left) = unsorted_rainflow(seq);

    // sort the cycles by peak position
    cycles.sort_by(|a, b| a.max.index.cmp(&b.max.index));

    (cycles, left)
}

/// find all the rainflow cycles in a sequence but do not sort them.
fn unsorted_rainflow<T: PartialOrd + Clone + std::fmt::Debug>(
    seq: &[T],
) -> (Vec<Cycle<T>>, Vec<T>) {
    // attach a label to the sequence
    let (mut cycles, mut left) = single_pass_rainflow(seq);
    let mut n = cycles.len();

    while n > 0 {
        let (mut new_cycles, new_left) = single_pass_rainflow(&left);
        n = new_cycles.len();
        cycles.append(&mut new_cycles);
        left = new_left;
    }

    (cycles, left)
}

///  Perform a single pass rainflow and extract all the cycles.
fn single_pass_rainflow<T: PartialOrd + Clone + std::fmt::Debug>(
    seq: &[T],
) -> (Vec<Cycle<T>>, Vec<T>) {
    let mut s = seq.iter();
    let mut cycles = Vec::new();
    let mut leftovers = Vec::new();

    let mut a = match s.next() {
        Some(x) => x,
        None => return (vec![], vec![]),
    };
    let mut b = match s.next() {
        Some(x) => x,
        None => return (vec![], vec![a.clone()]),
    };
    let mut c = match s.next() {
        Some(x) => x,
        None => return (vec![], vec![a.clone(), b.clone()]),
    };

    loop {
        leftovers.push(a.clone());
        let d = match s.next() {
            Some(x) => x,
            None => {
                // simple case of a three point sequence
                if a == c {
                    cycles.push(Cycle::new(b.clone(), c.clone()));
                } else {
                    leftovers.append(&mut vec![b.clone(), c.clone()]);
                }
                break;
            }
        };

        info!("between {:?} {:?} {:?} {:?}: {}", a,b,c,d, between_eq(a, b, d) && between_eq(a, c, d));

        if between_eq(a, b, d) && between_eq(a, c, d) {
            let cycle = Cycle::new(b.clone(), c.clone());

            cycles.push(cycle);
            a = d;
            b = match s.next() {
                Some(x) => x,
                None => {
                    leftovers.append(&mut vec![a.clone()]);
                    break;
                }
            };
            c = match s.next() {
                Some(x) => x,
                None => {
                    leftovers.append(&mut vec![a.clone(), b.clone()]);
                    break;
                }
            };
        } else {
            a = b;
            b = c;
            c = d;
        }
    }

    (cycles, leftovers)
}

/// Assign the cylces as simple tension cycles.
///
///  A tension cycle goes from a minimum to a maximum turning
///  point. Note: this routine needs some work to check the boundary
///  conditions at the start and end of the sequence.
pub fn tension<T: PartialOrd + Clone>(tp_sequence: &[T]) -> (Vec<Cycle<T>>, Vec<T>) {
    let mut tension_cycles = Vec::new();
    let mut unclosed = Vec::new();

    let start = if tp_sequence[0] < tp_sequence[1] {
        0
    } else {
        unclosed.push(tp_sequence[0].clone());
        1
    };

    let mut seq = start..(tp_sequence.len() - 1);
    while let Some(i) = seq.next() { 
        tension_cycles.push(Cycle::new(
            tp_sequence[i].clone(),
            tp_sequence[(i + 1)].clone(),
        ));
        seq.next();
    }

    (tension_cycles, unclosed)
}

/// Reorder the sequence to start at the maximum value.
pub fn reorder_sequence<T: PartialOrd + Clone>(seq: &[T]) -> Vec<T> {
    let mut start = 0;

    for (i, val) in seq.iter().enumerate() {
        if seq[start] < *val {
            start = i;
        }
    }
    let mut rseq = Vec::new();

    for i in start..seq.len() {
        rseq.push(seq[i].clone());
    }

    for i in 0..=start {
        rseq.push(seq[i].clone());
    }

    rseq
}

/// Print a summary of the list of cycles.
///
/// We cannot assume the cycles have come from rainflow counting so we
/// have a string containing the source of the cycles.
pub fn summarise_cycles(
    source: &str,
    cycles: &[Cycle<tag::Tag>],
    unclosed: &[tag::Tag],
    cycle_mods: &CycleModifiers,
) {
    println!(
        "
Cycle Summary
-------------
"
    );
    println!("Source: {}", source);
    println!("Cycle mods: {:#?}", cycle_mods);
    println!("Number of closed cycles: {}", cycles.len());
    println!("Number of unclosed turning points: {}", unclosed.len());

    if !cycles.is_empty() {
        let max = cycles
            .iter()
            .fold(cycles[0], |m, c| if c.range() > m.range() { *c } else { m });
        println!(
            "Largest range: {:.4e}, with valley of {} at line {}, peak of {} at line {}",
            max.range(),
            max.min.value,
            max.min.index,
            max.max.value,
            max.max.index
        );
        let min = cycles
            .iter()
            .fold(cycles[0], |m, c| if c.range() < m.range() { *c } else { m });
        println!(
            "Smallest range: {:.4e} with valley of {} at line {}, peak of {} at line {}",
            min.range(),
            min.min.value,
            min.min.index,
            min.max.value,
            min.max.index
        );
        println!(
            "Mean range: {:.4e}",
            cycles.iter().fold(0.0, |sum, x| sum + x.range()) / cycles.len() as f64
        );
        println!(
            "Maximum R: {:.4e}",
            cycles.iter().fold(cycles[0].r(), |m, p| p.r().max(m))
        );
        println!(
            "Minimum R: {:.4e}",
            cycles.iter().fold(cycles[0].r(), |m, p| p.r().min(m))
        );
        println!(
            "Mean R: {:.4e}",
            cycles.iter().fold(0.0, |sum, x| sum + x.r()) / cycles.len() as f64
        );
    }
}

/// Print a summary of the sequence.
pub fn summarise_sequence(source: &str, seq: &[tag::Tag], seq_mods: &SequenceModifiers) {
    println!(
        "
Sequence Summary
----------------
"
    );

    println!("Source: {}", source);

    // Write out all the modifications to the sequence.
    println!("Sequence mods: {:#?}", seq_mods);
    println!("Length: {}", seq.len());

    if !seq.is_empty() {
        let max = seq.iter()
            .fold(seq[0], |m, p| if p.value > m.value { *p } else { m });
        let max_count = seq.iter()
            .fold(0, |sum, s| sum + if (max.value - s.value).abs() < std::f64::EPSILON { 1 } else { 0 });
        println!("Maximum: {:.4e} at line {}", max.value, max.index);
        println!("Number of times maximum occurs: {}", max_count);
        let min = seq.iter()
            .fold(seq[0], |m, p| if p.value < m.value { *p } else { m });
        let min_count = seq.iter()
            .fold(0, |sum, s| sum + if (min.value - s.value).abs() < std::f64::EPSILON { 1 } else { 0 });
        println!("Minimum: {:.4e} at line {}", min.value, min.index);
        println!("Number of times minimum occurs: {}", min_count);
        println!(
            "Mean: {:.4e}",
            seq.iter().fold(0.0, |sum, x| sum + x.value) / seq.len() as f64
        );
        println!(
            "Number of points >= 0: {}",
            seq.iter()
                .fold(0, |count, x| if x.value >= 0.0 { count + 1 } else { count })
        );
        println!(
            "Number of points < 0: {}",
            seq.iter()
                .fold(0, |count, x| if x.value < 0.0 { count + 1 } else { count })
        );

        let tp = turning_points(seq);
        println!("Number of non turning points: {}", seq.len() - tp.len());
        print!("Sequence: ");
        print_start_end(seq, 5);
        println!()
    }
}

/// Print the n items from the start and end of a sequence in abbreviated form.
fn print_start_end(seq: &[tag::Tag], n: usize) {
    let len = seq.len();

    if len > 2 as usize * n {
        print!("[");
        for i in seq[0..n].iter() {
            print!("{:?} ", i.value);
        }

        print!("...");

        for i in seq[len - n..len].iter() {
            print!(" {:?}", i.value);
        }
        println!("]");
    } else {
        for i in seq.iter() {
            print!(" {:?}", i.value);
        }
    }
}

/// Extract the cycles from a sequence.
pub fn cycles_from_sequence(
    sequence: &[tag::Tag],
    cycle_method: &CycleMethod,
) -> (Vec<Cycle<tag::Tag>>, Vec<tag::Tag>) {
    let tp_reseq = turning_points(sequence);

    match *cycle_method {
        CycleMethod::Rainflow => rainflow(&tp_reseq),
        CycleMethod::Tension => tension(&tp_reseq),
    }
}

#[cfg(test)]
mod tests {
    use cycle::*;
    use tag;

    #[test]
    fn check_rainflow() {
        let s0 = vec![0.0, 1.0, 0.0];
        let ans0 = vec![Cycle::new(0.0, 1.0)];

        let s1 = vec![0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0];
        let ans1 = vec![Cycle::new(0.0, 1.0), Cycle::new(0.0, 2.0)];

        let s2 = vec![
            0.0, 10.0, 8.0, 5.0, 2.0, 1.0, 13.0, 8.0, 12.0, 8.0, 12.0, 6.0, 15.0, 6.0, 13.0, 0.0
        ];
        let ans2 = vec![
            Cycle::new(1.0, 10.0),
            Cycle::new(8.0, 12.0),
            Cycle::new(6.0, 13.0),
            Cycle::new(8.0, 12.0),
            Cycle::new(6.0, 13.0),
            Cycle::new(0.0, 15.0),
        ];

        let s3 = vec![0.0, 3.0, 2.0, 7.0, 4.0, 5.0, 4.0, 5.0, 0.0, 10.0, 0.0];
        let ans3 = vec![
            Cycle::new(3.0, 2.0),
            Cycle::new(5.0, 4.0),
            Cycle::new(10.0, 0.0),
            Cycle::new(5.0, 4.0),
            Cycle::new(7.0, 0.0),
        ];

        //        let ans3 = vec![];

        let s4 = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        let ans4 = vec![Cycle::new(1.0, 3.0)];

        let s5 = vec![0.0, -1.0, 1.0, 0.0];
        let ans5 = vec![];

        let s6 = vec![0.0, -1.0, 1.0, -1.0];
        let ans6 = vec![Cycle::new(-1.0, 1.0)];

        let s7 = vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let ans7 = vec![
            Cycle::new(0.0, 1.0),
            Cycle::new(0.0, 1.0),
            Cycle::new(0.0, 1.0),
        ];

        let s8 = vec![1.0, 3.0, 1.0];
        let ans8 = vec![Cycle::new(1.0, 3.0)];

        let s9 = vec![0.0, 1.0, 0.0, 1.0];
        let ans9 = vec![Cycle::new(0.0, 1.0)];

        let s10 = vec![0.0, 2.0, 1.0, 4.0, -2.0, 2.5, -2.5, 3.0, -3.0, 0.0];
        let ans10 = vec![
            Cycle::new(1.0, 2.0),
            Cycle::new(-2.0, 2.5),
            Cycle::new(-2.5, 3.0),
        ];

        // pathological sequence fails to produce any rainflow cycles
        let s11 = vec![0.0, 1.0, -2.0, 3.0, -4.0, 5.0, -3.0, 2.0, 1.0, 0.0];
        let ans11 = vec![];

        let all_seq = [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11];
        let all_ans = [
            ans0, ans1, ans2, ans3, ans4, ans5, ans6, ans7, ans8, ans9, ans10, ans11
        ];

        let tol = 1e-6;

        for (s, a) in all_seq.iter().zip(all_ans.iter()) {
            let tp = turning_points(s);
            // turn them into numbers not refs
            let tps = tp.iter().map(|&a| a).collect::<Vec<f64>>();

            let (rf, unclosed) = super::unsorted_rainflow(&tps);

            println!(
                "s {:?} tps {:?} rainflow {:?} unclosed {:?}",
                s, tps, rf, unclosed
            );
            assert!(rf.len() == a.len());
            for i in 0..rf.len() {
                assert!((rf[i].max - a[i].max).abs() < tol);
            }
        }
    }
    #[test]
    fn long_check_turning_points() {
        let mut x = Vec::new();
        let y = [0.0f64, 0.5, 1.0, 0.5];

        for _ in 0..100 {
            x.extend_from_slice(&y);
        }

        assert_eq!(turning_points(&x).len(), 201);
    }

    // #[test]
    // fn long_check_rainflow() {
    //     let mut x = Vec::new();
    //     let y = [0.0f64, 0.5, 0.25, 1.0, 0.5];

    //     for _ in 0..3 {
    //         x.extend_from_slice(&y);
    //     }

    //     let tp = turning_points(&x);
    //     let tps = tp.iter().cloned();
    //     let (rf, _unclosed) = super::unsorted_rainflow(&tps);
    //     println!("Unfiltered sequence {:?}", rf);

    //     // risefall and deadband filtering
    //     let risefall = rf.iter()
    //         .filter(|&cycle| cycle.range() > 0.5)
    //         .collect::<Vec<_>>();

    //     let deadband = rf.iter()
    //         .filter(|cycle| !cycle.within(0.5, 0.0))
    //         .collect::<Vec<_>>();

    //     let both = rf.iter()
    //         .filter(|&cycle| cycle.range() > 0.5)
    //         .filter(|cycle| !cycle.within(0.5, 0.0))
    //         .collect::<Vec<_>>();

    //     println!("risefall 0.5 sequence {:?}", risefall);
    //     println!("deadband [0.0, 0.5] sequence {:?}", deadband);
    //     println!("risefall and deadband sequence {:?}", both);

    // }

    #[test]
    fn check_tag_turning_points() {
        let s0 = tag::Tag::from(&[0.0, 1.0, 1.0, 4.0, 3.0, 5.0, 2.0, 0.0]);
        let ans0 = vec![
            tag::Tag::new(0.0, 0),
            tag::Tag::new(4.0, 3),
            tag::Tag::new(3.0, 4),
            tag::Tag::new(5.0, 5),
            tag::Tag::new(0.0, 7),
        ];

        let s1 = tag::Tag::from(&[1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0]);
        let ans1 = vec![
            tag::Tag::new(1.0, 0),
            tag::Tag::new(-1.0, 3),
            tag::Tag::new(1.0, 6),
        ];

        let s2 = tag::Tag::from(&[1.0, 0.0, 0.0, -1.0, -2.0, 0.0, 0.0, 1.0]);
        let ans2 = tag::Tag::from(&[1.0, -2.0, 1.0]);

        let s3 = tag::Tag::from(&[1.0, 1.5, 2.0, 2.5, 1.5, 1.0]);
        let ans3 = tag::Tag::from(&[1.0, 2.5, 1.0]);

        let s4 = tag::Tag::from(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        let ans4 = tag::Tag::from(&[1.0, 1.0]);

        let s5 = tag::Tag::from(&[1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0]);
        let ans5 = tag::Tag::from(&[1.0, 0.0, 1.0]);

        //    let s6 = tag::Tag::from(&[0.0]); // this should fail

        let seqs = [s0, s1, s2, s3, s4, s5];
        let ans = [ans0, ans1, ans2, ans3, ans4, ans5];

        // cargo test -- --nocapture
        for (sq, an) in seqs.iter().zip(ans.iter()) {
            let tp = turning_points(sq);
            println!("Seq: {:?}", sq);
            println!("Ans: {:?}", an);
            println!("Res: {:?}", tp);
            for (tpi, ai) in tp.iter().zip(an.iter()) {
                assert_eq!(*tpi, *ai);
            }
        }
    }
    #[test]
    fn check_cycle_order() {
        let c = Cycle::new(2.0, 1.0);
        assert!(c.max >= c.min);

        let c = Cycle::new(1.0, 2.0);
        assert!(c.max >= c.min);

        let c = Cycle::new(2.0, 2.0);
        assert!(c.max >= c.min);
    }

    #[test]
    fn check_simple_turning_points() {
        // the definitive turning point tests

        let s0 = vec![0.0, 1.0, 1.0, 4.0, 3.0, 5.0, 2.0, 0.0];
        let ans0 = vec![0.0, 4.0, 3.0, 5.0, 0.0];

        let s1 = vec![1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0];
        let ans1 = vec![1.0, -1.0, 1.0];

        let s2 = vec![1.0, 0.0, 0.0, -1.0, -2.0, 0.0, 0.0, 1.0];
        let ans2 = vec![1.0, -2.0, 1.0];

        let s3 = vec![1.0, 1.5, 2.0, 2.5, 1.5, 1.0];
        let ans3 = vec![1.0, 2.5, 1.0];

        let s4 = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let ans4 = vec![1.0, 1.0];

        let s5 = vec![1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0];
        let ans5 = vec![1.0, 0.0, 1.0];

        let s6 = vec![0.0];
        let ans6 = vec![0.0];

        let s7 = vec![0.0, 1.0];
        let ans7 = vec![0.0, 1.0];

        let s8 = vec![
            1.0, 1.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0
        ];
        let ans8 = vec![1.0, -1.0, 1.0];

        let all_s = [s0, s1, s2, s3, s4, s5, s6, s7, s8];
        let all_a = [ans0, ans1, ans2, ans3, ans4, ans5, ans6, ans7, ans8];

        // cargo test -- --nocapture
        for (s, a) in all_s.iter().zip(all_a.iter()) {
            let tp = turning_points(s);
            println!("seq: {:?} ans: {:?} tps:{:?}", s, a, tp);

            assert!(matching(&tp, a));
        }
    }

    fn matching<T: PartialEq>(a: &[T], b: &[T]) -> bool {
        a.iter().zip(b.iter()).filter(|&(a, b)| *a != *b).count() == 0
    }

    #[test]
    fn check_process_seq_mods() {
        let sequence = tag::Tag::from(&[0.0f64, 3.0, 2.0, 5.0, 3.0, 5.0, 4.0, 0.0]);

        let mods = SequenceModifiers {
            cap_max: Some(3.0),
            cap_min: None,
            remove_bigger: None,
            remove_smaller: None,
            cycles: false,
            reorder: false,
            turning_points: false,
            outfile: None,
        };

        let mod_sequence = process_seq_mods(&sequence, &mods);
        let correct_seq = tag::Tag::from(&[0.0f64, 3.0, 2.0, 3.0, 3.0, 3.0, 3.0, 0.0]);
        assert_eq!(correct_seq, mod_sequence);
    }
}
