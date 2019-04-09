use fatigue::{cycle, fracto, tag};
use options::{EasiOptions, OptimMethod, TerminatingOutput, Verbosity};
use nelder::Nelder;
use clap::{App, AppSettings, Arg};
use log::error;

/// Get the options from the command line.
pub fn get_options_clap(line: &str, options: &mut EasiOptions) {
    let process = App::new("easiGrow: A crack growth modelling program")
        .author("Paul White")
        .version(crate_version!())
        .about(include_str!("description.md"))
        .setting(AppSettings::AllowLeadingHyphen)
        
        .arg(Arg::with_name("scale")
             .short("s")
             .long("scale")
             .value_name("SCALE")
             .help("scale factor for load sequence (no default)")
             .takes_value(true))
        
        .arg(Arg::with_name("seq_infile")
             .short("q")
             .long("seq_infile")
             .value_name("FILE")
            .help("file of loading sequence for one block of loading")
             .takes_value(true))
        
        .arg(Arg::with_name("cycle_infile")
             .long("cycle_infile")
             .value_name("FILE")
             .help("read the cycles from an AFGROW input file")
             .takes_value(true))

        .arg(Arg::with_name("sequence")
             .long("sequence")
             .value_name("S1,S2,...")
             .help("explicitly set the sequence of loading (default: [0.0, 1.0, 0.0])")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("beta")
             .short("b")
             .long("beta")
             .value_name("NAME")
             .help("beta geometry model (default seft-newman84)")
             .takes_value(true))

        .arg(Arg::with_name("beta_outfile")
             .long("beta_outfile")
             .value_name("FILE")
             .help("write the beta model to a file")
             .takes_value(true))

        .arg(Arg::with_name("dadn")
             .short("d")
             .long("dadn")
             .value_name("NAME")
             .help("dadn material equation (default white:barter14_aa7050-t7451)")
             .takes_value(true))
                         
        .arg(Arg::with_name("params")
             .short("p")
             .long("parameters")
             .value_name("p1,p2,...,pM")
             .help("material parameters (default parameters for white:barter14_aa7050-t7451)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("cycle_method")
             .long("cycle_method")
             .help("method used to extract cycles from sequence (default rainflow)")
             .possible_values(&["rainflow", "tension"])
             .value_name("METHOD")
             .takes_value(true))

        .arg(Arg::with_name("astart")
             .short("a")
             .long("astart")
             .value_name("LENGTH")
             .help("initial crack sizes (default 10e-6 m). How many depends on beta factor, typically a,c is used.")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))
        
        .arg(Arg::with_name("forward")
             .long("forward")
             .value_name("DISTANCE")
             .help("distance forwards from crack origin to free edge in 'a' direction of growth (default: infinite)")
             .takes_value(true))
        
        .arg(Arg::with_name("sideways")
             .long("sideways")
             .value_name("DISTANCE")
             .help("distance sideways from crack origin to free edge in 'c' direction of growth (default: infinite)")
             .takes_value(true))
        
        .arg(Arg::with_name("radius")
             .long("radius")
             .value_name("R")
             .help("radius of hole or notch (default: infinite)")
             .takes_value(true))
        
    // Termination criteria 
    // only required if you want to terminate the crack calculation when this is exceeded
        .arg(Arg::with_name("aend")
             .short("e")
             .long("limit_a")
             .help("final termination crack sizes (default 1e-3 m)")
             .value_name("LENGTH")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("limit_k1c")
             .long("limit_k1c")
             .value_name("KIC")
             .help("fracture toughness K1C for termination (default 33 MPa sqrt(m))")
             .takes_value(true))
             
        .arg(Arg::with_name("limit_yield")
             .long("limit_yield")
             .value_name("STRESS")
             .help("yield or flow stress for plastic zone calculations and net section yield (default 450 MPa).")
             .takes_value(true))
        
        .arg(Arg::with_name("limit_block")
             .short("N")
             .long("limit_block")
             .value_name("N")
             .help("maximum number of blocks (default +1000)")
             .takes_value(true))

        .arg(Arg::with_name("youngs_modulus")
             .long("youngs")
             .value_name("MODULUS")
             .help("Young's modulus for compact tension coupon displacements (default 71000 MPa)")
             .takes_value(true))
        
    // Cycle counting
        .arg(Arg::with_name("seq_max")
             .long("seq_max")
             .value_name("MAX")
             .help("any value greater than this will be set to this value")
             .takes_value(true))

        .arg(Arg::with_name("seq_min")
             .long("seq_min")
             .value_name("MIN")
             .help("any value less than this will be set to this")
             .takes_value(true))

        .arg(Arg::with_name("seq_rem_big")
             .long("seq_rem_big")
             .value_name("BIG")
             .help("any values greater than this will be removed")
             .takes_value(true))

        .arg(Arg::with_name("seq_rem_small")
             .long("seq_rem_small")
             .value_name("SMALL")
             .help("any value less than this will be removed")
             .takes_value(true))

        .arg(Arg::with_name("seq_cycles")
             .long("seq_cyclemods")
             .help("only keep those turning points used in remaining cycles")
             .allow_hyphen_values(true))

        .arg(Arg::with_name("seq_outfile")
             .long("seq_outfile")
             .value_name("FILE")
             .help("write the modified sequence to a file")
             .takes_value(true))
                    
        .arg(Arg::with_name("cycle_max")
             .long("cycle_max")
             .value_name("MAX")
             .help("any range greater than this will be set to this value")
             .takes_value(true))

        .arg(Arg::with_name("cycle_min")
             .long("cycle_min")
             .value_name("MIN")
             .help("any range less than this will be set to this")
             .takes_value(true))

        .arg(Arg::with_name("cycle_rem_big")
             .long("cycle_rem_big")
             .value_name("DELTA")
             .help("remove cycles with a range bigger than this")
             .takes_value(true))

        .arg(Arg::with_name("cycle_rem_small")
             .long("cycle_rem_small")
             .value_name("DELTA")             
             .help("remove cycles with a range smaller than this")
             .takes_value(true))

        .arg(Arg::with_name("cycle_deadband")
             .long("cycle_deadband")
             .help("eliminates cycles with turning points within this band")
             .value_name("MIN,MAX")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("cycle_outfile")
             .long("cycle_outfile")
             .value_name("FILE")
             .help("write the modified sequence to a file")
             .takes_value(true))
                    
        .arg(Arg::with_name("seq_reorder")
             .short("r") 
             .long("seq_reorder")
             .help("re-order the sequence to close cycles (default: no reordering)"))

        .arg(Arg::with_name("seq_tp")
             .long("seq_tp")
             .help("remove non-turning points from the sequence (default: true)"))

        .arg(Arg::with_name("list")
             .short("l")
             .long("list")
             .help("list all available output variables, growth, beta and da/dN models"))

        .arg(Arg::with_name("output_every")
             .short("n")
             .long("output_every")
             .value_name("N")
             .help("output crack growth data every +N blocks or -N cycles (default +1)")
             .takes_value(true))

        .arg(Arg::with_name("output_lines")
             .long("output_lines")
             .value_name("l1,l2,...")
             .help("output growth at the given load lines l1,l2,l3... (default 1)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("output_vars")
             .short("o")
             .long("output_vars")
             .value_name("v1,v2,...")
             .help("comma separated list of variables to be output (default block,a)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))
                    
        .arg(Arg::with_name("verbose")
             .short("v")
             .long("verbose")
             .help("print more information"))

        .arg(Arg::with_name("summary")
             .long("summary")
             .help("print summary of sequence"))

        .arg(Arg::with_name("image_outfile")
             .long("image_outfile")
             .value_name("FILE")
             .help("generate a pseudo image of the fracture surface and write to FILE")
             .takes_value(true))
                    
        .arg(Arg::with_name("image_bar")
             .long("image_bar")
             .value_name("LENGTH")
             .help("size of scale bar for image (default 50e-6 m)")
            .takes_value(true))
        
        .arg(Arg::with_name("image_size")
             .long("image_size")
             .value_name("V,H")
             .help("size of image VxH pixels (default 300x8000)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("image_type")
             .long("image_type")
             .possible_values(&["Sem", "Optical"])
             .value_name("TYPE")
             .help("type of image output (default sem)")
            .takes_value(true))       

        .arg(Arg::with_name("opt_infile") 
             .long("opt_infile")
             .value_name("FILE")
             .help("optimises model parameters to match a set of target crack growth data.")
             .takes_value(true))
        
        .arg(Arg::with_name("opt_max") 
             .long("opt_max")
             .value_name("N")
             .help("maximum number of iterations to be used for optimisation (default 100)")
             .takes_value(true))

        .arg(Arg::with_name("opt_sweep")
             .long("opt_sweep")
             .value_name("f1,...,fM")
             .help("perform a brute force sweep over a range of scaling factors applied to M dadn parameters")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))
   
        .arg(Arg::with_name("opt_tol") 
             .long("opt_tol")
             .value_name("TOL")
             .help("terminate optimisation when the change in error function is less than this amount (default 1e-3)")
             .takes_value(true))

        .arg(Arg::with_name("opt_nelder") 
             .long("opt_nelder")
             .value_name("p1,...,p5")
             .help("nelder-mead parameters (default step: 0.1, alpha: 1.0, gamma: 2.0, rho: 0.5, sigma: 0.5)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("opt_method")
             .long("opt_method")
             .value_name("METHOD")
             .possible_values(&OptimMethod::variants())
             .help("method used for optimisation (default Levenberg)")
             .takes_value(true))
            
        .arg(Arg::with_name("crack_infile") 
             .short("c")
             .long("crack_infile")
             .value_name("FILE")
             .help("crack growth measurements that will be matched during optimisation")
             .takes_value(true))

        .arg(Arg::with_name("crack_weight") 
             .long("crack_weight")
             .value_name("WEIGHT")
             .help("weighting factor for matching crack growth curve (default 1.0)")
             .takes_value(true));

    // turn the commands into matches
    let matches = if line.is_empty() {
        // this will get them directly from the command line
        process.get_matches()
    } else {
        // get them from the string
        // add command word since it skips first argument
        process.get_matches_from(("easigrow ".to_string() + line.trim()).split(' '))
    };

    // basic crack growth options
    if let Ok(scale) = value_t!(matches, "scale", f64) {
        options.scale = scale;
    }
    if let Some(file) = matches.value_of("seq_infile") {
        options.seq_infile = file.to_string();
    }
    if let Some(file) = matches.value_of("cycle_infile") {
        options.cycle_infile = file.to_string();
    }
    if let Ok(sequence) = values_t!(matches, "sequence", f64) {
        options.sequence = tag::Tag::from(&sequence);
    }
    if let Some(beta) = matches.value_of("beta") {
        options.beta = beta.to_string();
    }
    if let Some(beta_outfile) = matches.value_of("beta_outfile") {
        options.beta_outfile = beta_outfile.to_string();
    }
    if let Some(dadn) = matches.value_of("dadn") {
        options.dadn = dadn.to_string();
    }
    if let Ok(params) = values_t!(matches, "params", f64) {
        options.params = params;
    }

    if let Some(cycle_method) = matches.value_of("cycle_method") {
        options.cycle_method = match cycle_method {
            "rainflow" => cycle::CycleMethod::Rainflow,
            "tension" => cycle::CycleMethod::Tension,
            _ => { error!("Error: unknown cycle counting method");
                   std::process::exit(2);
            },

        };
    }

    // crack and geometry options
    if let Ok(astart) = values_t!(matches, "astart", f64) {
        options.a = astart;
    }
    if let Ok(aend) = values_t!(matches, "aend", f64) {
        options.a_limit = aend;
    }
    if let Ok(forward) = value_t!(matches, "forward", f64) {
        options.component.forward = forward;
    }
    if let Ok(sideways) = value_t!(matches, "sideways", f64) {
        options.component.sideways = sideways;
    }
    if let Ok(radius) = value_t!(matches, "radius", f64) {
        options.component.radius = radius;
    }
    if let Ok(limit_k1c) = value_t!(matches, "limit_k1c", f64) {
        options.component.material.k1c = limit_k1c;
    }
    if let Ok(limit_yield) = value_t!(matches, "limit_yield", f64) {
        options.component.material.yield_stress = limit_yield;
    }
    if let Ok(limit_block) = value_t!(matches, "limit_block", f64) {
        options.block_limit = limit_block;
    }
    if let Ok(youngs_modulus) = value_t!(matches, "youngs_modulus", f64) {
        options.component.material.youngs_modulus = youngs_modulus;
    }

    // ssequence modifications
    if let Ok(max) = value_t!(matches, "seq_max", f64) {
        options.seq_mods.cap_max = Some(max);
    }
    if let Ok(min) = value_t!(matches, "seq_min", f64) {
        options.seq_mods.cap_min = Some(min);
    }
    if let Ok(big) = value_t!(matches, "seq_rem_big", f64) {
        options.seq_mods.remove_bigger = Some(big);
    }
    if let Ok(small) = value_t!(matches, "seq_rem_small", f64) {
        options.seq_mods.remove_smaller = Some(small);
    }
    if matches.is_present("seq_cycles") {
        options.seq_mods.cycles = true;
    }
    if let Some(output) = matches.value_of("seq_outfile") {
        options.seq_mods.outfile = Some(output.to_string());
    }
    if matches.is_present("seq_reorder") {
        options.seq_mods.reorder = true;
    }
    if matches.is_present("seq_tp") {
        options.seq_mods.turning_points = true;
    }

    // cycle modifications
    if let Ok(max) = value_t!(matches, "cycle_max", f64) {
        options.cycle_mods.cap_max = Some(max);
    }
    if let Ok(min) = value_t!(matches, "cycle_min", f64) {
        options.cycle_mods.cap_min = Some(min);
    }
    if let Ok(big) = value_t!(matches, "cycle_rem_big", f64) {
        options.cycle_mods.remove_bigger = Some(big);
    }
    if let Ok(small) = value_t!(matches, "cycle_rem_small", f64) {
        options.cycle_mods.remove_smaller = Some(small);
    }
    if let Ok(deadband) = values_t!(matches, "cycle_deadband", f64) {
        options.cycle_mods.remove_region = Some(cycle::Limit::from_vec(&deadband));
    }
    if let Some(output) = matches.value_of("cycle_outfile") {
        options.cycle_mods.outfile = Some(output.to_string());
    }

    // output options
    if let Ok(output_every) = value_t!(matches, "output_every", i32) {
        options.output_every = output_every;
    }
    if let Ok(output_vars) = values_t!(matches, "output_vars", String) {
        options.output_vars = output_vars;
    }
    if let Ok(output_lines) = values_t!(matches, "output_lines", usize) {
        options.output_lines = output_lines;
    }
    if matches.is_present("list") {
        options.output = TerminatingOutput::List;
    }
    if matches.is_present("summary") {
        options.output = TerminatingOutput::Summary;
    }
    if matches.is_present("verbose") {
        options.verbosity = Verbosity::Verbose;
        println!("Matches: {:?}", matches);
    }

    // Image generation options.
    if let Some(file) = matches.value_of("image_outfile") {
        options.image.file = file.to_string();
    }
    if let Ok(image_bar) = value_t!(matches, "image_bar", f64) {
        options.image.barlength = image_bar;
    }
    if matches.is_present("image_type") {
        options.image.image =
            value_t!(matches, "image_type", fracto::ImageType).unwrap_or_else(|e| e.exit())
    }
    if let Ok(image_size) = values_t!(matches, "image_size", u32) {
        options.image.xsize = image_size[0];
        options.image.ysize = image_size[1];
    }

    // Optimisation Options
    if let Some(file) = matches.value_of("opt_infile") {
        options.optimise.file = file.to_string();
    }
    if let Ok(opt_method) = value_t!(matches.value_of("opt_method"), OptimMethod) {
        options.optimise.method = opt_method;
    };
    if let Ok(opt_sweep) = values_t!(matches, "opt_sweep", f64) {
        options.optimise.sweep = opt_sweep;
    }
    if let Ok(opt_nelder) = values_t!(matches, "opt_nelder", f64) {
        options.optimise.nelder = Nelder::new(&opt_nelder);
    }
    if let Ok(opt_max) = value_t!(matches, "opt_max", usize) {
        options.optimise.maxiter = opt_max;
    }
    if let Ok(opt_tol) = value_t!(matches, "opt_tol", f64) {
        options.optimise.tol = opt_tol;
    }

    // Fracto options
    if let Some(file) = matches.value_of("crack_infile") {
        options.crack_infile = file.to_string();
    }
    if let Ok(crack_weight) = value_t!(matches, "crack_weight", f64) {
        options.crack_weight = crack_weight;
    }
}
