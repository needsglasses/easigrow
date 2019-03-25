//! Grow a fatigue crack from initial crack size until failure.

use std::f64::consts::FRAC_PI_2;
use COMMENT;
use std::process;
use std::f64::consts::PI;
use std::collections::BTreeSet;
use cycle::Cycle;
use plastic;
use dadn;
use material;
use tag;
use beta;

/// Data collected for each crack growth prediction cycle.
#[derive(Debug, Clone)]
pub struct History {
    /// Number of block.
    pub block: f64,
    /// Applied scaling stress for this cycle.
    pub stress: f64,
    /// Cycle information.
    pub cycle: Cycle<tag::Tag>,
    /// Stress intensity around the crack
    pub k: Vec<f64>,
    /// stress intensity range
    pub dk: Vec<f64>,
    /// beta values around the crack front
    pub beta: Vec<f64>,
    /// growth increment around the crack front
    pub da: Vec<f64>,
    pub crack: CrackState,
}

#[derive(Debug, Clone)]
pub struct CrackState {
    /// Crack length at each point around the crack front
    pub a: Vec<f64>,
    /// distance of monotonic zone from crack tip
    pub mono_zone_extent: plastic::ZoneWidth,
    /// distance of cyclic zone from crack tip
    pub cyclic_zone_extent: plastic::ZoneWidth,
}

impl CrackState {
    pub fn new(a: Vec<f64>) -> CrackState {
        CrackState {
            a,
            mono_zone_extent: plastic::ZoneWidth {
                plane_stress: 0.0,
                plane_strain: 0.0,
            },
            cyclic_zone_extent: plastic::ZoneWidth {
                plane_stress: 0.0,
                plane_strain: 0.0,
            },
        }
    }
}

#[derive(Debug, Clone)]
/// describes the geometry and material containing the crack
pub struct Component {
    pub forward: f64,
    pub sideways: f64,
    pub radius: f64,
    pub material: material::Properties,
}

// All the data necessary to run a fatigue calculation.
pub struct FatigueTest {
    pub history: History,
    pub component: Component,
    pub scale: f64,
    pub cycles: Vec<Cycle<tag::Tag>>,
    pub a_limit: Vec<f64>,
    pub block_limit: f64,
    pub next_cycle: usize,
    pub dadn: Box<dadn::DaDn>,
    pub beta: Box<beta::Beta>,
    pub output_vars: Vec<String>,
}

// Create an interator that can step through a fatigue test struct to
// generate the history for each cycle.  needs a 'history' and a set
// of 'options'.
impl Iterator for FatigueTest {
    type Item = History;

    fn next(&mut self) -> Option<History> {
        // increment the next cycle
        self.history.block =
            self.history.block.floor() + (self.next_cycle as f64 + 1.0) / self.cycles.len() as f64;
        self.next_cycle = (self.next_cycle + 1) % self.cycles.len();
        let tagged_cycle = self.cycles[self.next_cycle];

        let history = self.history.clone();
        // grow the crack for one cycle
        self.history = history.grow_crack(
            &tagged_cycle,
            self.scale,
            &self.dadn,
            &mut self.beta,
            &self.component,
        );

        // Check for terminating conditions.
        let reached_limit = reached_limit(
            self.history.block,
            self.block_limit,
            &self.history.crack.a,
            &self.a_limit,
        );
        let component_failed = component_failed(
            &self.history.crack.a,
            self.history.stress * self.history.cycle.max.value,
            self.history.k[0],
            &self.component,
            &self.beta,
        );

        if reached_limit.failure || component_failed.failure {
            // Print the final line at failure.
            display_history_line(&self.history, &self.output_vars, &self.component);
            let messages = reached_limit.messages + &component_failed.messages;
            println!("{}", messages);
            None
        } else {
            Some(self.history.clone())
        }
    }
}

/// Performs a crack growth calculation for a single cycle.
impl History {
    pub fn grow_crack<T: dadn::DaDn + ?Sized, U: beta::Beta + ?Sized>(
        self,
        cycle: &Cycle<tag::Tag>,
        scale: f64,
        eqn: &Box<T>,
        beta: &mut Box<U>,
        component: &Component,
    ) -> History {
        // grow the crack cycle by cycle
        let Cycle {
            max: tag::Tag {
                value: smax,
                ..
            },
            min: tag::Tag {
                value: smin,
                ..
            },
        } = *cycle;

        let r = smin / smax;
        if smax < smin {
            println!(
                "Program Error: smax {} is less than smin {}. This should never happen.",
                smax, smin
            );
            process::exit(1);
        }

        let c = self.crack.a[self.crack.a.len() - 1];
        let a_on_c = self.crack.a[0] / c;
        let a_on_d = self.crack.a[0] / component.forward;
        let c_on_b = c / component.sideways;
        let a_on_r = self.crack.a[0] / component.radius;
        let phis = vec![0.0, FRAC_PI_2];

        let betas = beta.beta(a_on_d, a_on_c, c_on_b, a_on_r, &phis);

        // values around the crack front
        let mut da_all: Vec<f64> = Vec::new();
        let mut k_all: Vec<f64> = Vec::new();
        let mut dk_all: Vec<f64> = Vec::new();

        //        println!("Scale is {}", scale);
        // calculate the growth for each 'a' around the crack
        // front. We assume the beta is calculated at each point
        // around the crack front but we know that for some beta
        // functions this is not possible.
        let mut a_all = Vec::with_capacity(self.crack.a.len());
        for (beta, a) in betas.iter().zip(self.crack.a.iter()) {
            let k_on_stress = scale * beta * (PI * *a).sqrt();
            let kmin = smin * k_on_stress;
            let kmax = smax * k_on_stress;
            let dk = kmax - kmin;

            k_all.push(kmax);
            dk_all.push(dk);

            // grow the crack an increment
            let da = eqn.dadn(kmin, kmax, dadn::CrackState{ a: *a });
            if da.is_nan() {
                println!(
                    "Error: the dadn calculation has returned an NAN beta {} cycle_dk {} r {}",
                    beta, dk, r
                );
                process::exit(1);
            }
            a_all.push(a + da);
            da_all.push(da);
        }

        // size of the plastic zone
        let kmax = k_all.iter().fold(k_all[0], |f, x| x.max(f));
        let dkmax = dk_all.iter().fold(dk_all[0], |f, x| x.max(f));

        let mono_zone_extent = plastic::zone_size(kmax, component.material.yield_stress);
        // Simply double the yield stress to get the cyclic stress yield stress
        let cyclic_zone_extent = plastic::zone_size(dkmax, 2. * component.material.yield_stress);

        History {
            block: self.block, // part_block,
            da: da_all.clone(),
            k: k_all.clone(),
            dk: dk_all.clone(),
            cycle: cycle.clone(),
            stress: scale,
            beta: betas,
            crack: CrackState {
                a: a_all,
                mono_zone_extent,
                cyclic_zone_extent,
            },
        }
    }
}

pub fn display_history_header(output: &[String]) {
    // write out the headers
    if !output.is_empty() {
        for out in output.iter() {
            print!("{:>12} ", out);
        }
        println!();
    }
}

/// Test to see if a history line should be output.
pub fn output_cycle_history(
    his: &History,
    every: i32,
    output_lines: &BTreeSet<usize>,
    cycle_no: usize,
) -> bool {
    // if every is positive write out every nth block, otherwise if
    // every is negative write out every nth cycle
    let frequency = (every > 0 && his.block as i32 % every == 0)
        || (every < 0 && cycle_no % -every as usize == 0);
    // println!("freq {} every {}, lines {:?}, block {}, cycle_no {}, max {} min {}", frequency, every, output_lines, his.block, cycle_no, his.cycle.max.index,his.cycle.min.index);

    // output only if the cycle constains the specific sequence line
    if !output_lines.is_empty() && every > 0 {
        frequency
            && (output_lines.contains(&his.cycle.max.index)
                || output_lines.contains(&his.cycle.min.index))
    } else {
        frequency
    }
}

// print a line of the history data
pub fn display_history_line(his: &History, output: &[String], component: &Component) {
    let a = his.crack.a[0];
    let c = his.crack.a[his.crack.a.len() - 1];

    for out in output {
        match out.trim() {
            "line" => print!("{:12} ", his.cycle.max.index),
            "block" => print!("{:12.4} ", his.block),
            "a/c" => print!("{:12.4} ", a / c),
            "a/d" => print!("{:12.4} ", a / component.forward),
            "c/b" => print!("{:12.4} ", c / component.sideways),
            "k" => print!("{:12.4} ", his.k[0]),
            "dk" => print!("{:12.4} ", his.dk[0]),
            "r" => print!("{:12.4} ", his.cycle.min.value / his.cycle.max.value),
            "beta_a" => print!("{:12.4e} ", his.beta[0]),
            "beta_c" => print!("{:12.4e} ", his.beta[his.beta.len() - 1]),
            "a" => print!("{:12.6e} ", a),
            "c" => print!("{:12.6e} ", c),
            "da" => print!("{:12.4e} ", his.da[0]),
            "dc" => print!("{:12.4e} ", his.da[his.da.len() - 1]),
            "mono" => print!("{:12.4e} ", his.crack.mono_zone_extent.plane_strain),
            "cyclic" => print!("{:12.4e} ", his.crack.cyclic_zone_extent.plane_strain),
            "a/mono" => print!("{:12.2} ", a / his.crack.mono_zone_extent.plane_strain),
            "a/cyclic" => print!("{:12.4} ", a / his.crack.cyclic_zone_extent.plane_strain),
            "mono/da" => print!(
                "{:12.4}",
                his.crack.mono_zone_extent.plane_strain / his.da[0]
            ),
            "cyclic/da" => print!(
                "{:12.4}",
                his.crack.cyclic_zone_extent.plane_strain / his.da[0]
            ),
            "peak" => print!("{:12.6} ", his.stress * his.cycle.max.value),
            "valley" => print!("{:12.6} ", his.stress * his.cycle.min.value),
            ref opt => {
                println!(
                    "Error: Unknown output option (use the --list option for a complete list): {}",
                    opt
                );
                process::exit(1);
            }
        };
    }
    println!();
}

/// This structure provides a way of remembering the failure messages
/// so that they can be written out at the end.
pub struct FailureResult {
    /// Type of failure.
    pub failure: bool,
    /// Failure Message.
    pub messages: String,
}

/// check if we have reached a pre-defined crack growth limit.
pub fn reached_limit(
    part_block: f64,
    block_limit: f64,
    a: &[f64],
    a_limit: &[f64],
) -> FailureResult {
    let mut message = "".to_string();
    let mut terminate = false;

    if part_block >= block_limit {
        message += &format!(
            "{}Run stopped because hit block limit {}\n",
            COMMENT, block_limit
        );
        terminate = true;
    }

    if !a_limit.is_empty() && a.iter().zip(a_limit).any(|(a, e)| a >= e) {
        message += &format!(
            "{}Failure Event: a{:?} >= a_limit{:?}\n",
            COMMENT, a, a_limit
        );
        terminate = true;
    }

    FailureResult {
        failure: terminate,
        messages: message,
    }
}

/// Check whether the crack has exceeded any failure criteria for the component.
pub fn component_failed(
    a: &[f64],
    smax: f64,
    kmax: f64,
    component: &Component,
    _beta: &Box<beta::Beta>,
) -> FailureResult {

    let mut terminate = false;
    let mut message = "".to_string();

    // clippy says to do it this way, but its a bit inconsistent with the multiple failure checks
    // Check whether we have satisfied any termination criteria.
    // let mut terminate = if component.forward > 0.0 && a[0] > component.forward {
    //     message += &format!(
    //         "{}Failure Event: a[{}] > depth[{}]\n",
    //         COMMENT, a[0], component.forward
    //     );
    //     true } else { false };

    if component.forward > 0.0 && a[0] > component.forward {
        message += &format!(
            "{}Failure Event: a[{}] > depth[{}]\n",
            COMMENT, a[0], component.forward
        );
        terminate = true;
    }

    if component.material.k1c > 0.0 && kmax > component.material.k1c {
        message += &format!(
            "{}Failure Event: k[{}] > k1c[{}]\n",
            COMMENT, kmax, component.material.k1c
        );
        terminate = true;
    }

    // The net stress for the component will depend on the shape of a crack.
    // We assume it is a corner crack (but the worst case will be an internal crack).
    let approx_crack_area = (PI / 4.0) * a.first().unwrap() * a.last().unwrap();
    let component_area = component.sideways * component.forward;
    let applied_stress = smax * component_area / (component_area - approx_crack_area);

    if component.material.yield_stress > 0.0 && applied_stress > component.material.yield_stress {
        message += &format!(
            "{}Note: Assuming a corner crack to check the net-section yield stress criterion.\n",
            COMMENT
        );
        message += &format!(
            "{}approx crack area {}, component area {}\n",
            COMMENT, approx_crack_area, component_area
        );
        message += &format!(
            "{}Failure Event: stress[{}] > yield[{}]\n",
            COMMENT, applied_stress, component.material.yield_stress
        );
        terminate = true;
    }

    FailureResult {
        failure: terminate,
        messages: message,
    }
}
