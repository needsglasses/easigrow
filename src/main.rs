/// easiGro
///
/// by Paul White (Nov 2014--2017)
/// written in rust (www.rust-lang.org)
///
/// A program to match crack growth predictions to measurements.
///
/// The program calculates fatigue crack growth rates and finds the
/// optimum parameters of a crack growth model to match predictions
/// with measurements.
///
/// **easiGrow** is a standalone program but most of the calculations
/// are done through calls to the associated **fatigue** library which
/// is included. The main program is for doing anything that
/// explicitly uses the command line flags inlcuding the optimisation
/// module. These flages are used to build the **EasiOptions** data
/// structure which is then used to generate the crack growth
/// history. The optimisation generates a crack growth curve which it
/// compares with a fractography file. It finds the error between
/// these measurements and tries to minimise the sum errors through
/// minimisation routines.
///
/// Currently, none of the models has a memory effect, so it is ok to
/// just start growing the crack from an iniital crack size that is
/// smaller than the initial fracto data. The struct `grow::CrackState`
/// also contains parameters that are passed along with the applied
/// loading _kmin_ and _kmax_, so any memory variables should be added to
/// this struct and will be availabe to be used by the _da/dn_ equation.
/// The simplest memory effect that is included in the `CrackState`
/// data is the plastic zone size, but there are no dadn equations
/// currently using this. The memory effect does not appear to be
/// strong in AA7050 material.
///
/// Think of the program flow as
///
/// 1. Read in data
/// 2. Filter the sequence (turning point, rainflow, risefall, deadband etc.) and convert to cycles
/// 3. Filter the list of cycles
/// 4. If required, optimise any parameters
/// 5. Perform a crack growth calculation
/// 6. Write out requested output

// #![cfg_attr(feature="clippy", feature(plugin))]
// #![cfg_attr(feature="clippy", plugin(clippy))]
use std::f64::consts::FRAC_PI_2;
use std::process;
use std::collections::BTreeSet;
use fatigue::{beta, cycle, dadn, fracto, grow, io, material, table, tag};
use options_clap::get_options_clap;
use options::{OptimMethod, TerminatingOutput};
use fatigue::dadn::DaDn;
use fatigue::COMMENT;
use std::error::Error;
use std::fs::File;
use std::path::Path;
use log::error;
use std::io::Write;

extern crate log;
extern crate env_logger;

#[macro_use]
extern crate clap;
extern crate fatigue;

//mod options_argparse;
mod list;
mod optimise;
mod sweep;
mod factors;
mod options;
mod options_clap;
mod nelder;
mod numbers;
mod vector;

#[cfg(feature = "GSL")]
mod optimise_gsl;

fn main() {
    env_logger::init();
    
    // get all the data
    let materials = material::get_all_dadns();

    let mut options = options::get_default_options();
    get_options_clap("", &mut options);
    println!("{}easiGrow: version {}", COMMENT, crate_version!());
    println!("{}", COMMENT);
    if options.verbosity == options::Verbosity::Verbose {
        println!("{}Options: ", COMMENT);
        println!("{}", options);
    }
    
    options::read_all_files(&mut options);

    // process all the modifications to the sequence
    options.sequence = cycle::process_seq_mods(&options.sequence, &options.seq_mods);

    // Get the cycles from either the external sequence file, command line or the cycle file.
    if options.cycle_infile != "" && options.seq_infile != "" {
        error!("Error: you have specified a sequence file '{}' as well as a cycle file '{}'. Specify only one.",
                 options.seq_infile, options.cycle_infile);
        std::process::exit(2)
    }

    let unclosed = if options.cycle_infile == "" {
        let (cycles, left) = cycle::cycles_from_sequence(&options.sequence, &options.cycle_method);
        options.cycles = cycles;
        left
    } else {
        Vec::new()
    };    

    // process all the modifications to the cycles
    options.cycles = cycle::process_cycle_mods(&options.cycles, &options.cycle_mods);

    // Only keep those cycles that remain after filtering the cycles
    // and mark the turning points associated with those cycles. This
    // section is only for writing out the modified sequence, since
    // the filtered cycles are all that is used for crack growth.
    if options.seq_mods.cycles {
        let mut keep = vec![false; options.sequence.len()];
        for cycle in &options.cycles {
            keep[cycle.max.index] = true;
            keep[cycle.min.index] = true;
        }
        options.sequence.retain(|s| keep[s.index])
    }

    // Any request for file or info output will result in program
    // termination. This policy is to reduce the complexity for the
    // user as to what the program does.

    // Write out the sequence file.
    if let Some(outfile) = options.seq_mods.outfile {
        io::write_sequence(&outfile, &options.sequence);
        std::process::exit(0);
    }

    // Write out the cycles file.
    if let Some(outfile) = options.cycle_mods.outfile {
        io::write_cycles(&outfile, &options.cycles);
        std::process::exit(0);
    }

    // write out the beta by converting to a beta table. This can be
    // then read back in using the file: option for beta selection.
    if options.beta_outfile != "" {
        let beta = beta::get_beta_fn(&options.beta, &options.component);
        let table_beta = beta.as_table();

        // need to write to file
        let path = Path::new(&options.beta_outfile);
        let display = path.display();

        let mut file = match File::create(&path) {
            // The `description` method of `io::Error` returns a string that
            // describes the error
            Err(why) => {
                error!(
                    "Error: could not create the file '{}': {}.",
                    display,
                    Error::description(&why)
                );
                std::process::exit(1)
            }
            Ok(file) => file,
        };
        let _ = write!(file, "{}", table_beta);
        std::process::exit(0);
    }

    // write out summary information of the sequence
    match options.output {
        TerminatingOutput::Summary => {
            let seq_source = if options.seq_infile != "" {
                options.seq_infile
            } else {
                // This is a little vague as the sequence could be either
                // the default sequence or overwritten with a supplied sequence.
                String::from("(Used specified sequence)")
            };
            cycle::summarise_sequence(&seq_source, &options.sequence, &options.seq_mods);

            let cycle_source = if options.cycle_infile != "" {
                options.cycle_infile
            } else {
                format!(
                    "(Obtained from sequence using '{:?}' method)",
                    options.cycle_method
                )
            };
            cycle::summarise_cycles(
                &cycle_source,
                &options.cycles,
                &unclosed,
                &options.cycle_mods,
            );

            std::process::exit(0)
        }

        // write out extended list of options and methods
        TerminatingOutput::List => {
            list::print_list();
            std::process::exit(0);
        }

        _ => (),
    }

    // get the correct material parameters for the dadn equation or
    // from the command line. If the params are not given, then get the
    // dadn material constants from the internal database.
    let mut params = options.params.clone();
    if params.is_empty() {
        // extract out the appropriate material parameters from a file
        params = if options.dadn.starts_with("file:") {
            let filename = options.dadn.trim_start_matches("file:");
            println!(
                "{}No parameters given, using the dk values in the dadn file {}",
                COMMENT, filename
            );
            let table = table::Table::read_file(filename, true);
            // collapse down the dks and use these as the parameters for optimising
            table.variables()
        // or from the internal database.
        } else {
            println!(
                "{}No parameters given, obtaining from material library for {}",
                COMMENT, options.dadn
            );
            match materials.iter().find(|m| options.dadn.starts_with(m.name)) {
                Some(m) => m.eqn.variables(),
                None => {
                    error!("Error: Unknown dadn model {}", options.dadn);
                    process::exit(1);
                }
            }
        }
    };

    // Optimise the parameters to match the predicted crack growth
    // rates with the associated measured crack growth rates.
    if options.optimise.file != "" {
        // optimisation scaling factors
        options.params = params.clone();
        println!(
            "{}Now starting the optimisation with params {:?} ...",
            COMMENT, options.params
        );

        let mut factors = vec![1.0; params.len()]; // non-dimensionalised factors used for optimisation

        optimise_error(&options, &mut factors);
        println!("{}...finished the optimisation. ", COMMENT);
        println!("{}The normalised factors are {:?}", COMMENT, factors);

        // Rescale the parameters to include the optimised factors
        params = options
            .params
            .iter()
            .zip(factors)
            .map(|(p, f)| p * f)
            .collect::<Vec<f64>>();

        println!("{}The scaled optimised factors are: {:?}", COMMENT, params);

        if options.scale == 0.0 {
            std::process::exit(0); // not an error if we have performed an optimisation
        }
    }

    // Grow the crack
    let history_all = generate_crack_history(&options, &params);

    // Lastly, now that we've grown the crack, check if we need to
    // generate and write out a pseudo image.
    if options.image.file != "" {
        println!("Making a pseudo image...");
        if options.image.file.ends_with(".svg") {
            fracto::write_svg_pseudo_image(&history_all, &options.image, &options.image.file);
            println!("image written to file {}", &options.image.file);
        } else {
            error!("Error: Currently easigo can only generate svg. Please use a '.svg' suffix");
        }
    }
}

#[cfg(not(feature = "GSL"))]
fn optimise_error(options: &options::EasiOptions, mut factors: &mut [f64]) {
    match options.optimise.method {
        OptimMethod::Sweep => sweep::sweep(options, &mut factors),
        OptimMethod::Nelder => optimise::nelder_match_crack(options, &mut factors),
        OptimMethod::All => {
            sweep::sweep(options, &mut factors);
            optimise::nelder_match_crack(options, &mut factors)
        }
    };
}

#[cfg(feature = "GSL")]
fn optimise_error(options: &options::EasiOptions, mut factors: &mut [f64]) {
    match options.optimise.method {
        OptimMethod::Sweep => sweep::sweep(options, &mut factors),
        OptimMethod::Nelder => optimise::nelder_match_crack(&options, &mut factors),
        OptimMethod::Levenberg => optimise_gsl::gsl_match_crack(&options, &mut factors),
        OptimMethod::All => {
            sweep::sweep(options, &mut factors);
            optimise::nelder_match_crack(options, &mut factors);
            optimise_gsl::gsl_match_crack(options, &mut factors)
        }
    };
}

// Finally grow the crack with the current parameters which may have been optimised.

// We exit here if the scale has not been set. Otherwise we
// would go through and do a default calculation which confuses
// people if they just want to start the program to see how to get
// help.
fn generate_crack_history(options: &options::EasiOptions, params: &[f64]) -> Vec<grow::History> {
    let dadn_eqn = dadn::make_model(&options.dadn, &params, String::from("unknown"));
    println!("{}da/dN equation: {}", COMMENT, dadn_eqn);

    let beta = beta::get_beta_fn(&options.beta, &options.component);

    if options.scale == 0.0 {
        error!(
            "Error: The sequence scale factor is 0. You need to set the scale factor
           (i.e. load or stress level) in order to perform a crack growth calculation.
           Try\n easigrow --help"
        );
        std::process::exit(1);
    }

    if options.cycles.is_empty() {
        println!("Error: There are no closed cycles in sequence. Perhaps try the re-order sequence option -r");
        std::process::exit(1);
    }

    // We define the initial state. If any memory effect is to be
    // included in the crack growth model, the meory should be in this
    // data structure.
    let init_crack = grow::CrackState::new(options.a.clone());

    let mut history_all = Vec::new();
    grow::display_history_header(&options.output_vars);

    // Non-dimensional ratios for beta factor
    let c = options.a[options.a.len() - 1];
    let a_on_c = options.a[0] / c;
    let a_on_d = options.a[0] / options.component.forward;
    let c_on_b = c / options.component.sideways;
    let a_on_r = options.a[0] / options.component.radius;
    let phis = vec![0.0, FRAC_PI_2];

    // Initialise the history
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
        beta: beta.beta(a_on_d, a_on_c, c_on_b, a_on_r, &phis),
        da: vec![0.0, 0.0],
        crack: init_crack,
    };

    grow::display_history_line(&init_history, &options.output_vars, &options.component);
    let component = grow::FatigueTest {
        history: init_history,
        component: options.component.clone(),
        scale: options.scale,
        cycles: options.cycles.clone(),
        a_limit: options.a_limit.clone(),
        block_limit: options.block_limit,
        next_cycle: 0,
        dadn: dadn_eqn,
        beta,
        output_vars: options.output_vars.clone(),
    };

    // make a hash set of the lines that are required for output
    let mut output_lines: BTreeSet<usize> = options.output_lines.iter().cloned().collect();

    // if there are no lines in the output then put in the line for the first cycle
    if options
        .cycles
        .iter()
        .filter(|c| output_lines.contains(&c.max.index) || output_lines.contains(&c.min.index))
        .count() == 0
    {
        println!("output_lines {:?}", output_lines);
        println!(
            "
Warning: There are no sequence lines in the cycle list and so there
         will be no crack growth output. Consider closing up cycles
         with re-order to use all sequence lines or include specific
         sequence lines that are in the cycle. Meanwhile, the output will
         be for the squence line in the first cycle at line {}.",
            options.cycles[0].max.index
        );
        output_lines.insert(options.cycles[0].max.index);
    }

    // Start the crack growth. This loop steps through each cycle
    // repeating the cycles until a terminating condition stops the
    // growth and ends the for loop.
    for (cycle_no, history) in component.enumerate() {
        if grow::output_cycle_history(&history, options.output_every, &output_lines, cycle_no) {
            grow::display_history_line(&history, &options.output_vars, &options.component);
        }
        // Only keep the history if we are producing a fracto image.
        if options.image.file != "" {
            history_all.push(history);
        }
    }
    history_all
}
