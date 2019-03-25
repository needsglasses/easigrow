//! This is done by a call to the minpack Levenberg-Marquardt routine
//! implemented in GSL that minimises a vector of differences of the
//! growth rates. The growth rates are calculated from the target
//! crack growth files by taking the differences between successive
//! crack measurements in the file.

use std::process;
use options;
use optimise;

extern crate rgsl;

#[derive(Clone)]
struct Data {
    options: Vec<options::EasiOptions>,
    params: Vec<f64>,
}

// Optimise the parameters of the crack growth model to minimise the
// error in matching crack growth data. Factors are the scaling
// factor applied to the equation parameters that are to be optimised.
pub fn gsl_match_crack(main_options: &options::EasiOptions, factors: &mut [f64]) -> f64 {
    let all_options = optimise::get_all_options(main_options, &main_options.optimise.file);

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

    let data = Data {
        options: all_options,
        params: main_options.params.clone(),
    };
    println!("GSL Params {:?}", data.params);
    // the optimised variables are multiplied by the parameters
    let x = rgsl::types::vector::VectorF64::from_slice(factors).unwrap();

    rgsl::RngType::env_setup();
    let t = rgsl::MultiFitFdfSolverType::lmsder();

    let mut f = rgsl::MultiFitFunctionFdf::new(m, n);

    // this is a copy of the code in gsl examples
    // let predict_f = clone!(data => move |x, f| {
    //     println!("In the new lambda: {:?}", x);
    //     gsl_prediction_error(&x, &mut f, &*data.borrow())
    // });

    let predict_f = move |x, mut f| {
        println!("In the lambda: {:?}", x);
        gsl_prediction_error(&x, &mut f, &data)
    };

    f.f = Some(Box::new(predict_f));

    // let mut f = rgsl::MultiFitFunctionFdf {
    //     f: gsl_prediction_error,
    //     df: None,
    //     fdf: None,
    //     n: m,
    //     p: n,
    //     params: &mut data
    // };

    let mut iter = 0;
    let mut s = rgsl::MultiFitFdfSolver::new(&t, m, n).unwrap();

    s.set(&mut f, &x);

    //    print_state(iter as u64, &mut s);
    let mut status = rgsl::Value::Continue;

    loop {
        println!("looping");
        if iter >= main_options.optimise.maxiter {
            // don't warn for single iteration passes
            if (main_options.verbosity == options::Verbosity::Verbose)
                && main_options.optimise.maxiter > 0
            {
                println!("Warning: we have exceeded the specified maximum iteration limit. Iter = {} > {} {:?}", 
                         iter, main_options.optimise.maxiter, status);
            }
            break;
        }

        iter += 1;

        println!("before s.x: {:?}", s.x());
        status = s.iterate();
        println!("after s.x: {:?}", s.x());

        if status != rgsl::Value::Success {
            break;
        }

        if main_options.verbosity == options::Verbosity::Verbose {
            println!("Checking convergence change in s.dx {:?}", s.dx());
            println!("Checking convergence change in s.x {:?}", s.x());
        }

        status = rgsl::multifit::test_delta(&s.dx(), &s.x(), 1e-10, 1e-10);
        if status != rgsl::Value::Continue {
            break;
        }
    }

    if main_options.verbosity == options::Verbosity::Verbose {
        println!("status = {}", rgsl::error::str_error(status));
        println!("Finished the optimisation");
    }

    for i in 0..x.len() {
        factors[i as usize] = s.x().get(i);
    }

    rgsl::blas::level1::dnrm2(&s.f())
}

// Wrapper around prediction error to make it compatible with gsl
fn gsl_prediction_error(x: &rgsl::VectorF64, f: &mut rgsl::VectorF64, data: &Data) -> rgsl::Value {
    let mut params = Vec::new();

    // x is the vector of nondimensionalised parameters
    for (i, p) in data.params.iter().enumerate() {
        params.push(p * x.get(i).max(0.0).min(2.0));
    }
    println!("gsl_prediction_error: {:?}", params);
    // calculate the error term for each measurement difference
    let errors = optimise::prediction_error(&params, &data.options);

    // add all the individal error measurements into a global error vector
    let mut j = 0;
    for growth_error in errors {
        for g in growth_error {
            f.set(j, g);
            j += 1;
        }
    }

    rgsl::Value::Success
}
