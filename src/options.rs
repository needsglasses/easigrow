/// These are the data stuctures for command line options as well as
/// the initial default values.

use std;
use std::f64;
use std::string::String;
use nelder::Nelder;
use fatigue::{cycle, fracto, grow, io, material, tag};
use fatigue::COMMENT;

#[cfg(not(feature = "GSL"))]
arg_enum!{
    #[derive(Debug, Clone)]
    pub enum OptimMethod {
        Sweep,
        Nelder,
        All
    }
}

#[cfg(feature = "GSL")]
arg_enum!{
    #[derive(Debug, Clone)]
    pub enum OptimMethod {
        Sweep,
        Nelder,
        Levenberg,
        All
    }
}

arg_enum!{
    #[derive(Debug, Clone)]
    pub enum CycleMethod {
        Rainflow,
        Tension
    }
}

#[derive(Debug, Clone)]
pub struct Optimise {
    pub file: String,
    pub method: OptimMethod,
    pub maxiter: usize,
    pub tol: f64,
    pub nelder: Nelder,
    pub sweep: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Verbosity {
    Verbose,
    Terse,
}
#[derive(Debug, Clone, PartialEq)]
pub enum TerminatingOutput {
    List,
    Summary,
    None,
}

/// Option data for performing a crack growth calculation.
#[derive(Debug, Clone)]
pub struct EasiOptions {
    /// Initial starting size, last value in vector is considered as c the crack width.
    pub a: Vec<f64>,
    /// Final crack depth.
    pub a_limit: Vec<f64>,
    /// Maximum number of cycles to run for.
    pub block_limit: f64,
    /// Type of output
    pub output: TerminatingOutput,
    /// Level of verbosity
    pub verbosity: Verbosity,
    /// Name of inbuilt dadn data to use.
    pub dadn: String,
    /// Dadn parameters to be used in optimisation or any crack calculation.
    pub params: Vec<f64>,
    /// Closure parameters
    pub closure: Vec<f64>,
    /// Parameters to be written out.
    pub output_vars: Vec<String>,
    /// Output 'every' block.
    pub output_every: i32,
    /// Sequence lines to write out
    pub output_lines: Vec<usize>,
    /// Name of damage model to use.
    pub cycle_method: cycle::CycleMethod,
    /// Scale the sequence by this factor, typically stress
    pub scale: f64,
    /// Beta model to use.
    pub beta: String,
    pub beta_outfile: String,
    /// Information describing the shape of the component.
    pub component: grow::Component,
    /// Sequence to be used.
    pub sequence: Vec<tag::Tag>,
    /// Simplified version of the measured crack file.
    pub fracto: Vec<io::Measurement>,
    /// Data from crackfile.
    pub cycles: Vec<cycle::Cycle<tag::Tag>>,
    /// Info for optimisation.
    pub optimise: Optimise,
    /// Data for generating a reconstructed image.
    pub image: fracto::ImageData,
    /// Name of crack growth data file used for target in optimisation.
    pub crack_infile: String,
    // Weighting factor for crack errors
    pub crack_weight: f64,
    /// Sequence information.
    pub seq_mods: cycle::SequenceModifiers,
    pub cycle_mods: cycle::CycleModifiers,
    /// Filename of sequence.
    pub seq_infile: String,
    /// Name of cycle file for inputting data.
    pub cycle_infile: String,
}


impl std::fmt::Display for EasiOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let _e = writeln!(f, "{}a: {:?}", COMMENT, self.a);
        let _e = writeln!(f, "{}a_limit: {:?}", COMMENT, self.a_limit);
        let _e = writeln!(f, "{}block_limit: {:?}", COMMENT, self.block_limit);
        let _e = writeln!(f, "{}params: {:?}", COMMENT, self.params);
        let _e = writeln!(f, "{}output_vars: {:?}", COMMENT, self.output_vars);
        let _e = writeln!(f, "{}output_every: {:?}", COMMENT, self.output_every);
        let _e = writeln!(f, "{}output_lines: {:?}", COMMENT, self.output_lines);
        let _e = writeln!(f, "{}scale: {:?}", COMMENT, self.scale);
        let _e = writeln!(f, "{}beta: {:?}", COMMENT, self.beta);
        let _e = writeln!(f, "{}component: {:?}", COMMENT, self.component);
        let _e = writeln!(f, "{}seq_infile: {:?}", COMMENT, self.seq_infile);
        let _e = writeln!(f, "{}seq_mods: {:?}", COMMENT, self.seq_mods);
        let _e = writeln!(f, "{}cycle_method: {:?}", COMMENT, self.cycle_method);
        let _e = writeln!(f, "{}cycle_infile: {:?}", COMMENT, self.cycle_infile);
        let _e = writeln!(f, "{}cycle_mods: {:?}", COMMENT, self.cycle_mods);

        if self.optimise.file != "" {
            let _e = writeln!(f, "{}optimise: {:?}", COMMENT, self.optimise);
            let _e = writeln!(f, "{}crack_infile: {:?}", COMMENT, self.crack_infile);
            let _e = writeln!(f, "{}crack_weight: {:?}", COMMENT, self.crack_weight);
        }

        if self.image.file != "" {
            let _e = writeln!(f, "{}fracto: {:?}", COMMENT, self.fracto);
            let _e = writeln!(f, "{}image: {:?}", COMMENT, self.image);
        }
        write!(f, "{}dadn: {:?}", COMMENT, self.dadn)
    }
}

/// read in each sequence and the measured crack growth file
///
/// This stuff only needs to be performed once such as populating the
/// values of tables read from files.
pub fn read_all_files(options: &mut EasiOptions) {
    if !options.seq_infile.is_empty() {
        options.sequence = io::read_sequence(&options.seq_infile);
    }

    if !options.cycle_infile.is_empty() {
        options.cycles = io::read_afgrow_cycles(&options.cycle_infile);
    }

    if !options.crack_infile.is_empty() {
        options.fracto = io::read_fracto_file(&options.crack_infile, options.sequence.len());
    }
}

pub fn get_default_options() -> EasiOptions {
    // default options for easigro
    EasiOptions {
        // crack growth formulation
        a: vec![10_e-6f64, 10e-6],
        dadn: "white:barter14-aa7050t7451".to_string(),
        params: vec![],
        closure: vec![0.3, 0.5],
        beta: "seft-newman84".to_string(),
        beta_outfile: "".to_string(),
        component: grow::Component {
            sideways: f64::INFINITY,
            forward: f64::INFINITY,
            radius: f64::INFINITY,
            material: material::Properties {
                yield_stress: 450.0,
                k1c: 33.0,
                youngs_modulus: 71e3,
            },
        },

        // termination criteria
        a_limit: vec![1e-3, 1e-3], // (m)
        block_limit: 1000.0,

        // sequence info
        scale: 0.0,                                 // MPa
        seq_infile: "".to_string(),                 // name of sequence file
        sequence: tag::Tag::from(&[0.0, 1.0, 0.0]), // Constant amplitude cycle
        cycles: vec![],                             // cycles from
        cycle_infile: "".to_string(),
        cycle_method: cycle::CycleMethod::Rainflow,
        seq_mods: cycle::SequenceModifiers {
            cap_max: None,
            cap_min: None,
            remove_bigger: None,
            remove_smaller: None,
            cycles: false,
            reorder: false,
            turning_points: false,
            outfile: None,
        },
        cycle_mods: cycle::CycleModifiers {
            cap_max: None,
            cap_min: None,
            remove_bigger: None,
            remove_smaller: None,
            remove_region: None,
            outfile: None,
        },

        // crack comparison
        crack_infile: "".to_string(), // name of fracto file
        crack_weight: 1.0,
        fracto: vec![],

        // solution type
        optimise: Optimise {
            file: "".to_string(),
            maxiter: 100,
            method: OptimMethod::Nelder,
            tol: 1e-4,
            nelder: Nelder::default(),
            sweep: vec![0.8, 1.0, 1.25],
        },

        image: fracto::ImageData {
            file: "".to_string(),
            barlength: 50e-6,
            xsize: 300,
            ysize: 8000,
            image: fracto::ImageType::Sem,
        },

        // output options
        output_vars: vec!["block".to_string(), "a".to_string()],
        output_every: 1, // output 'every' block
        output_lines: vec![1],
        output: TerminatingOutput::None,
        verbosity: Verbosity::Terse,
    }
}
