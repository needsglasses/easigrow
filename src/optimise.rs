//! This module optimises the model to best match crack growth measurements.
//!

use std::error::Error;
use std::f64;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process;
use nelder;
use options;
use options_clap::get_options_clap;
use fatigue::{beta, cycle, dadn, grow, plastic, tag};
use self::rayon::prelude::*;
use fatigue::io::Measurement;
use std::collections::BTreeMap;

extern crate rayon;

/// Match crack growth measurements using Nelder-Mead algorthm.
pub fn nelder_match_crack(main_options: &options::EasiOptions, factors: &mut [f64]) -> f64 {
    let all_options = get_all_options(main_options, &main_options.optimise.file);

    // number of variables to be optimised
    let n = factors.len(); // number of optimise parameters

    // Calculate the total number of points to be matched across all
    // match files. This may be excessive and could be reduced to one
    // value for each measured crack curve.
    let m: usize = all_options.iter().fold(0, |sum, option| {
        let nrun = 1.0f64;
        sum + (option.fracto.len() as f64 - nrun) as usize
    });

    if main_options.verbosity == options::Verbosity::Verbose {
        println!("match_crack: n (variables) {} m (targets) {}", n, m);
    }

    if m < n {
        println!(
            "Error: Insufficient number of match points {} to optimise {} variables",
            m, n
        );
        process::exit(1);
    }

    let sum_error_fn = |x: &[f64]| sum_prediction_error(x, &all_options);

    println!("Initial x: {:?}, f(x): {}", factors, sum_error_fn(factors));

    nelder::nelder_mead(
        &sum_error_fn,
        factors,
        &main_options.optimise.nelder,
        main_options.optimise.tol,
        main_options.optimise.maxiter,
    )
}

/// Read all the option lines from the optimisation file.
pub fn get_all_options(
    default_options: &options::EasiOptions,
    file: &str,
) -> Vec<options::EasiOptions> {
    let path = Path::new(file);
    let handle = match File::open(&path) {
        // The `description` method of `io::Error` returns a string that
        // describes the error
        Err(why) => {
            println!(
                "Error: could not open file '{}': {}.",
                path.display(),
                why.description()
            );
            process::exit(1)
        }
        Ok(file) => file,
    };

    let optimfile = BufReader::new(handle);

    // create a vector of the process options for each optimisation
    let mut all_options: Vec<options::EasiOptions> = Vec::new();

    // read in each line of the optimisation file and store as option arguments
    for line in optimfile.lines() {
        let mut options = default_options.clone();

        match line {
            Ok(line) => {
                if default_options.verbosity == options::Verbosity::Verbose {
                    println!("matchfile: {:?}", line);
                }
                get_options_clap(&line, &mut options);
            }
            Err(e) => {
                println!(
                    "Error: problem in reading line from the optimisation file: {}.",
                    e
                );
                process::exit(1);
            }
        }

        // make sure we are using the main dadn model in every sub model
        options.dadn = default_options.dadn.clone();
        println!("crackfile: {}", options.crack_infile);
        options::read_all_files(&mut options);

        // process all the modifications to the cycles
        options.sequence = cycle::process_seq_mods(&options.sequence, &options.seq_mods);

        if options.cycle_infile != "" && options.seq_infile != "" {
            println!("Error: you have specified a sequence file '{}' as well as a cycle file '{}'. Specify only one.",
                     options.seq_infile, options.cycle_infile);
            process::exit(2)
        }

        let (cycles, _unclosed) =
            cycle::cycles_from_sequence(&options.sequence, &options.cycle_method);
        // process all the modifications to the cycles
        options.cycles = cycle::process_cycle_mods(&cycles, &options.cycle_mods);

        if options.crack_infile.is_empty() {
            println!(
                "Error: The --crack_infile option is missing from the list of
            parameters in the optimisation file. The crackfile option
            is necessary for optimisation because it supplies the target
            crack growth curve."
            );
            process::exit(1);
        }

        all_options.push(options);
    }
    all_options
}

/// Calculates the total error between predictions and measurements.
pub fn sum_prediction_error(x: &[f64], options: &[options::EasiOptions]) -> f64 {
    let mut params = Vec::new();
    // x is a nondimensionalised parameter
    for (i, p) in options[0].params.iter().enumerate() {
        params.push(p * x[i].max(0.0)); // restrict variable from going negative
    }

    let errors = prediction_error(&params, options);
    let mut sum = 0.0;

    for (e, o) in errors.iter().zip(options) {
        let mut err = 0.0;
        for f in e {
            err += f * f;
        }
        sum += err * o.crack_weight
    }
    sum.sqrt()
}

/// Error function is the difference between predicted and measured crack growth rates.
///
/// This is typically the objective function to be minimised.
pub fn prediction_error(params: &[f64], options: &[options::EasiOptions]) -> Vec<Vec<f64>> {
    let mut errors: Vec<Vec<f64>> = Vec::with_capacity(options.len());

    // calculate the errors
    println!("Using the material parameters {:?}", params);

    options
        .par_iter()
        .map(|option| {
            // find the lines in the measurements that have been used
            // so that we only check the corresponding line numbers in the history
            let mut lines_used = BTreeMap::new();
            //            println!("fracto: {:?}", option.fracto);
            for meas in &option.fracto {
                lines_used.entry(meas.line).or_insert(meas.line);
            }
            lines_used
                .entry(option.cycles[0].max.index)
                .or_insert(option.cycles[0].max.index);
            println!("optimising by comparing lines: {:?}", lines_used);
            let init_crack = grow::CrackState {
                a: option.a.clone(),
                mono_zone_extent: plastic::ZoneWidth {
                    plane_stress: 0.0,
                    
                    plane_strain: 0.0,
                },
                cyclic_zone_extent: plastic::ZoneWidth {
                    plane_stress: 0.0,
                    plane_strain: 0.0,
                },
            };
            let init_history = grow::History {
                block: 0.0,
                stress: 0.0,
                cycle: cycle::Cycle {
                    max: tag::Tag {
                        value: 0.0,
                        index: 0,
                    },
                    min: tag::Tag {
                        value: 0.0,
                        index: 0,
                    },
                },
                k: vec![0.0, 0.0],
                dk: vec![0.0, 0.0],
                beta: vec![0.0, 0.0],
                da: vec![0.0, 0.0],
                crack: init_crack,
            };

            let dadn_eqn = dadn::make_model(&option.dadn, &params, String::from("unknown"));
            let component = grow::FatigueTest {
                history: init_history,
                component: option.component.clone(),
                scale: option.scale,
                cycles: option.cycles.clone(),
                a_limit: option.a_limit.clone(),
                block_limit: option.block_limit,
                next_cycle: 0,
                dadn: dadn_eqn,
                beta: beta::get_beta_fn(&option.beta, &option.component),
                output_vars: option.output_vars.clone(),
            };

            // Only keep in history those sequence lines that are used in the fracto measurements
            // this is an approximation since there may be nearer points that are not included.
            let history = component
                .filter(|h| {
                    // println!("his: {:?}", h);
                    lines_used.contains_key(&h.cycle.max.index)
                        || lines_used.contains_key(&h.cycle.min.index)
                })
                .collect::<Vec<_>>();

            compare_growth(&option.fracto, &history, &option.verbosity)
        })
        .collect_into_vec(&mut errors);

    let mut total_error = Vec::new();
    for (option, prediction_error) in options.iter().zip(&errors) {
        let rms_error = rms(prediction_error);
        println!("{}: {}", option.crack_infile, rms_error);
        total_error.push(rms_error);
    }

    println!("Total error: {}", rms(&total_error));

    errors
}

/// Calculate the Root-Mean-Square value of a sequence.
pub fn rms(x: &[f64]) -> f64 {
    x.iter().fold(0.0f64, |acc, &a| acc + a.powi(2)).sqrt()
}

/// Calculate the errors between predicted and measured data.
fn compare_growth(
    measured: &[Measurement],
    history: &[grow::History],
    verbosity: &options::Verbosity,
) -> Vec<f64> {
    let mut errors = Vec::with_capacity(measured.len());

    // println!("comparing the growth");
    // println!("history {}", history.len());

    // for h in history {
    //     println!("{} {:12.6e}", h.block, h.crack.a[0]);
    // }

    // println!("Measured {}", measured.len());
    // for m in measured {
    //     println!("{} {:12.6e}", m.block, m.a);
    // }

    // Find the predicted crack growth rates at the measured crack
    // size by looking for the corresponding history that match the
    // two adjacent fracto measurements.
    for i in 0..(measured.len() - 1) {
        // Skip the comparison if crack length is zero implying a
        // discontinuous measurement.
        if (measured[i].a == 0.0) || (measured[i + 1].a == 0.0) {
            continue;
        }

        // Check through the history for each measurement to find a
        // predicted crack size that is bigger than the measured size.
        let h_end = match history.iter().skip(1).position(|his| {
            his.crack.a[0] >= measured[i + 1].a && if measured[i + 1].line == 0 {
                true
            } else {
                measured[i + 1].line == his.cycle.max.index
            }
        }) {
            Some(h_end) => h_end + 1,
            // can't find one so use the last point
            None => history.len() - 1,
        };

        // Similarly, start searching backwards through the history to
        // find the latest predicted size that is smaller than the measured size.
        let h_start = match history.iter().rev().skip(1)
        // don't let the start be the last point otherwise we will get a zero block length
            .position(|his| his.crack.a[0] < measured[i].a 
                      && his.crack.a[0] < history[h_end].crack.a[0]
                      && if measured[i+1].line == 0 {
                          true
                      } else {
                          measured[i+1].line == his.cycle.max.index
                      }) {
            Some(h_start) => history.len() - h_start - 2,
            None => 0,
        };

        // check if we do not grow the crack from a small enough size
        if history.is_empty() {
            panic!("Zero length crack growth history");
        }

        if history[0].crack.a[0] - history[0].da[0] > measured[i].a {
            println!(
                "Warning: predicted starting crack size {} is greater than measured data {}",
                history[0].crack.a[0] - history[0].da[0],
                measured[i].a
            );
        }

        let history_dadb = (history[h_end].crack.a[0] - history[h_start].crack.a[0])
            / (history[h_end].block - history[h_start].block);
        let measured_dadb =
            (measured[i + 1].a - measured[i].a) / (measured[i + 1].block - measured[i].block);

        let error = if (history.is_empty()) || (history_dadb == 0.0) {
            println!(
                "Warning: history block size is zero h_start {} a_start {} h_end {} a_end {}.",
                h_start, history[h_start].crack.a[0], h_end, history[h_end].crack.a[0]
            );
            -100.0
        } else {
            // The choice of error term is difficult. Here we use the
            // difference between logs, which says that the order of
            // magnitude of the accuracy is important with equal
            // weighting being given to being out by an order of
            // magnitude at high crack growth rates as with small
            // crack growth rates. In the absence of any knowledge on
            // the number of cycles of each this seems like the best
            // error function.
            history_dadb.ln() - measured_dadb.ln()
        };

        if *verbosity == options::Verbosity::Verbose {
            println!("New Measurement comparison summary:");
            println!(
                "  History: start {:6.2} at a={}, end {:6.2} at a={}",
                history[h_start].block,
                history[h_start].crack.a[0],
                history[h_end].block,
                history[h_end].crack.a[0]
            );
            println!(
                "  Measurement: start {} at a={}, end {} at a={}",
                measured[i].block,
                measured[i].a,
                measured[i + 1].block,
                measured[i + 1].a
            );
            println!(
                "  giving average growth rates: history {:8e}, measured {:8e}, error {:8e}",
                history_dadb, measured_dadb, error
            );
        }

        if h_end == h_start {
            println!("Programming error: Things have gone wrong in compare_growth() with hstart {} == hend {}",
                     h_start, h_end);
        }

        errors.push(error);
    }

    errors
}

#[cfg(test)]
mod tests {
    use super::compare_growth;
    use fatigue::io::Measurement;
    use fatigue::{cycle, grow, plastic, tag};
    use options;

    #[test]
    fn test_compare_growth() {
        let measured = vec![
            Measurement {
                line: 0,
                block: 20.0,
                a: 0.11,
            },
            Measurement {
                line: 0,
                block: 30.0,
                a: 0.20,
            },
            Measurement {
                line: 0,
                block: 40.0,
                a: 0.31,
            },
        ];

        let history = vec![
            fake_history(110.0, 0.05),
            fake_history(120.0, 0.15),
            fake_history(130.0, 0.25),
            fake_history(140.0, 0.35),
            fake_history(150.0, 0.45),
        ];

        let errors = compare_growth(&measured, &history, &options::Verbosity::Verbose);

        let mut e = errors.into_iter();

        assert!((e.next().unwrap() - 0.1054).abs() < 1e-3);
        assert!((e.next().unwrap() - -0.09531).abs() < 1e-3);
    }

    fn fake_history(block: f64, a: f64) -> grow::History {
        grow::History {
            block,
            k: vec![10.0, 20.0],
            dk: vec![5.0, 10.0],
            stress: 1.0,
            cycle: cycle::Cycle {
                min: tag::Tag::new(5.0, 0),
                max: tag::Tag::new(10.0, 0),
            },
            beta: vec![0.1, 0.2],
            da: vec![0.1, 0.3],
            crack: grow::CrackState {
                a: vec![a, 0.2],
                mono_zone_extent: plastic::ZoneWidth {
                    plane_stress: 0.0,
                    plane_strain: 0.0,
                },
                cyclic_zone_extent: plastic::ZoneWidth {
                    plane_stress: 0.0,
                    plane_strain: 0.0,
                },
            },
        }
    }

}
