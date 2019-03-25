use options;
use factors;
use optimise;
use numbers::NonNan;
use std::{process, f64};
use self::rayon::prelude::*;

extern crate rayon;

/// Perform a brute force sweep over all permutations of the parameters.
///
/// Calculates the error function for each permutation of the
/// factors.  This function can be used as a brute force optimisation to search
/// a space for the best fit.  Or, it could be used for a one off error
/// function evaluation for all of the crack files that are to be matched.
pub fn sweep(options: &options::EasiOptions, params: &mut [f64]) -> f64 {
    println!("Performing a brute force sweep. This may take some time...");
    let sweep_factors = factors::permutations(&options.optimise.sweep, params.len());
    println!("Sweep: there are {} combinations", sweep_factors.len());
    println!("Sweep: parameters {:?}", options.params);

    if options.verbosity == options::Verbosity::Verbose {
        for factor in &sweep_factors {
            println!("Sweep factors {:?}", factor);
        }
    }
    if options.optimise.file == "" {
        println!(
            "Error: you need an optimisation file to perform a sweep in order to calculate errors"
        );
        process::exit(1);
    }

    let all_options = optimise::get_all_options(options, &options.optimise.file);
    let mut results = Vec::with_capacity(sweep_factors.len());
    sweep_factors
        .par_iter()
        .map(|normalised_factors| {
            let scaled_factors = normalised_factors
                .iter()
                .zip(options.params.iter())
                .map(|(&f, &p)| f * p)
                .collect::<Vec<_>>();
            let single_error = optimise::prediction_error(&scaled_factors, &all_options);
            let total_single_error = single_error
                .iter()
                .map(|x| optimise::rms(x).powi(2))
                .sum::<f64>()
                .sqrt();
            (normalised_factors, total_single_error)
        })
        .collect_into(&mut results);

    for &(f, r) in &results {
        println!("{:?} {:?}", f, r);
    }

    // print out the results for the smallest error
    let (best_normalised_factors, smallest_error) = results
        .into_iter()
        .min_by_key(|&(_factor, error)| NonNan::new(error.abs()))
        .unwrap();

    println!("Best result from sweep:");
    println!("    Total Error: {:?}", smallest_error);
    println!("    Normalised factors: {:?}", best_normalised_factors);

    // copy back the best parameters
    for i in 0..params.len() {
        params[i] = best_normalised_factors[i];
    }

    smallest_error
}
