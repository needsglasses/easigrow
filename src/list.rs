use fatigue::{beta, grow, material};

static HIGHLIGHTS: &'static str = "
* **Sequence filtering** Performs sequence reordering,
   turning-point, dead-band, rise-fall filtering and rain-flow
   counting.  The filtered sequences may be written to a file.

* **Inbuilt data** Comes with a selection of beta factors,
   material data and crack growth models.

* **Calculated parameters** Calculates additional parameters for
   characterising the state of the crack tip so that better
   crack growth equations can be developed based on the most
   applicable parameters for a material.

* **Optimisation** Optimises the crack growth model parameters to
   minimise the difference between predicted and measured crack
   growth. The measured data need not be for the entire history
   i.e., one or more fractographic measurements of the width of
   a block.

* **Image generation** Generates a pseudo fractographic image of the fracture
   surface to see how easy it is to identify the individual
   blocks. 
";

static UNITS: &'static str = "The internal da/dN data are all in units
for stress intensity of (MPa m^0.5) and growth in (m). Most beta
equations use an applied far field stress which is in units of
(MPa). Thus the scaling factor in units of (MPa) will generally
convert the sequence into the correct units. However, the
compact-tension beta factor compact_tada73 uses applied load not
stress and is in units of load of (MN) and coupon dimensions are in
(m). The width and depth will need to be set for the compact-tension
beta function otherwise **easiGrow** will assume an infinite plate and
the crack will not grow.";

/// Prints out lists of data. Sort of an extended help.
pub fn print_list() {
    // List of the valid names for producing output
    let output = [
        ("block", "block number"),
        ("line", "line number"),
        ("a/c", "crack aspect ratio"),
        ("a/d", "cracked fraction of forward distance"),
        ("c/b", "cracked fraction of sideways distance"),
        ("k", "current cycle K"),
        ("dk", "current cycle dK"),
        ("r", "current cycle R"),
        ("beta_a", "beta factor at a"),
        ("beta_c", "beta factor at c"),
        ("a", "forward distance of crack from origin"),
        ("c", "sideways distance of crack from origin"),
        ("da", "crack growth increment at a"),
        ("dc", "crack growth increment at c"),
        ("mono", "largest monotonic plastic zone "),
        ("cyclic", "largest cyclic plastic zone "),
        (
            "a/mono",
            "ratio of crack length to monotonic plastic zone size",
        ),
        (
            "a/cyclic",
            "ratio of crack length to cyclic plastic zone size",
        ),
        (
            "mono/da",
            "ratio of current cyclic plastic zone to current da",
        ),
        (
            "cyclic/da",
            "ratio of current cyclic plastic zone to current da",
        ),
        ("peak", "scaled peak stress of current cycle"),
        ("valley", "scaled valley stress of current cycle"),
    ];

    let formats = [
        (
            "Crack file",
            "The crack growth file is in the following format:

<line> <block> <a>
...

or 

<block> <a>
...

or

<a>
...
                   
Blank lines in the file indicate non-contiguous measurements. If
<line> or <block> are missing the program will assume the readings are
one block apart with each block measured at line 0. Use the same
format for the entire file. Where <line> represents the corresponding
line no. (starting at 0) of the sequence file, and <block> is the
no. of the block at that crack depth. Strictly speaking, the actual
block numbers are not used by easigrow with only the difference between
the block numbers in contiguous measurements used. Easigrow only
matches the average crack growth rate using:

   rate = log(growth between measurements) / log(no. of blocks between measurements).
",
        ),
        (
            "Optimise file",
            "Each line in the optimise file is a list of easigrow command lines
(without the 'easigrow' command) with each line containing the easigrow
options that will best reproduce the crack growth calculation for the
associated crack growth curve that it is trying to match. Note: Only
the material model specified by the main command line that invokes the
optimisation will be used for all crack predictions, since those will
be the parameters that are optimised. Any other material
specifications will be ignored.

The format of the optimisation file is:

<easigrow option> ... --crack <FILE1>
<easigrow option> ... --crack <FILE2>
...
",
        ),
        (
            "Beta file",
            "All lines beginning with a # are treated as a comment and ignored. The
format of the beta file is

# Comment describing the contents of the file
a/d beta 
...
",
        ),
        (
            "Dadn file",
            "All lines beginning with a # are treated as a comment and ignored. The
format of the file is:

# Comment describing the contents of the file
r1 r2 ....
dadn1 deltaK1_r1 deltaK1_r2 ....
dadn2 deltaK2_r1 deltaK2_r2 ....
...
",
        ),
    ];

    let biblio = [
        ["[Newman79]", "J. C. Newman , Jr. and I. S. Raju
           Analyses of surface cracks in finite plates under tension or bending loads
           NASA Technical Paper 1578
           December 1979"],
        
        ["[Newman81]", " J. C. Newman Jr. and I. S. Raju
           Stress intensity factor equations for cracks
           in three-dimensional finite bodies, 
           NASA Technical Memorandum 83299, 1981 p 1--49."],

        ["[Newman81]", " J. C. Newman Jr. and I. S. Raju
           Stress-intensity factor equations for
           cracks in three-dimensional finite bodies
           subjected to tension and bending loads, 
           NASA Technical Memorandum 85739, 1984."],

        ["[Anderson05]", "T. L. Anderson
           Fracture Mechanics - Fundamentals and Applications
           Taylor and Francis 3rd Edition 2005"],
        
        ["[Tada73]", "H. Tada, P. C. Paris and G. R. Irwin
           The Stress Analysis of Cracks Handbook
           1973"],
        
        ["[Murakami87]", "Y. Murakami
          Stress Intensity Factors Handbook. Vol 2
          Pergamon Press, Oxford, , 1987"],

        ["[Murakami87a]", "Yukitaka Murakami and Hideto Tsuru
          Stress Intensity factor equations for a semi-elliptical surface crack in a shaft under bending
          1986"],

        ["[Schwarmann86]", "L. Schwarmann
          Material Data of High-Strength Aluminium Alloys for Durability Evaluation of Structures
          Aluminium-Verlag 1986
          Note: The data from this report has been converted from mm/cycle to m/cyclic by factoring cf by 1e3."],

        ["[Fedderson66]", "
         Taken from Damage Tolerant Design handbook from AFGROW documentation."],

        ["[Kujawski01]", "Daniel Kujawski, 
          A fatigue crack driving force parameter with load ratio effects
          International Journal of Fatigue, Vol 23, S239-S246, 2001"],

        ["[Walker70]", "K. Walker
          The effect of stress ratio during crack propagation and fatigue for {2024-T3} and {7075-T6} aluminum
          Effects of Environment and Complex Load History for Fatigue Life, 
          American Society for Testing and Materials,Special Technical Publication 462, 1970"],

        ["[Jones12]", "Jones, R., Molent, L. & Walker, K.  
          Fatigue crack growth in a diverse
          range of materials, International Journal of Fatigue Vol. 40,pages 43--50, 2012"],
        
        ["[Hartman70]", "A. Hartman and J. Schijve
          The effects of environment and load frequency on the
          crack propagation law for macro fatigue crack growth in aluminum alloys,
          Engineering Fracture Mechanics, Vol. 1(4), 1970"],

        ["[Shin04]", "C.S. Shin and C. Q. CAI
          Experimental and finite element analyses on stress intensity
          factors of an elliptical surface crack in a circular shaft under
          tension and bending,
          International Journal of Fracture 129: 239–264, 2004."],

        ["[Forman05]", "R. G. Forman, V. Shivakumar, J. W. Cardinal , L. C. Williams and P. C. McKeighan 
                        Fatigue Crack Growth Database for Damage Tolerance Analysis, 
                        DOT/FAA/AR-05/15, 2005."],
    ];

    // Set up a new counter to automatically label the section headers
    let mut header = Counter::new();

    header.section("Program Highlights");
    print!("{}", HIGHLIGHTS);

    header.section("Units");
    print!("{}", UNITS);

    header.section("Output Parameters");

    for &(abbrev, descrip) in &output {
        println!("{:20} {}", abbrev, descrip);
    }

    header.section("Beta Models");
    // List of the beta equations that are implemented
    let betas = beta::get_all_betas(&grow::Component {
        forward: 0.0,
        sideways: 0.0,
        radius: 0.0,
        material: material::Properties {
            yield_stress: 0.0,
            k1c: 0.0,
            youngs_modulus: 0.0,
        },
    });
    for beta in &betas {
        println!(
            "{:20} {:30} {} {}",
            beta.name,
            beta.args.to_string(),
            beta.summary,
            beta.cite
        );
    }

    header.section("Cycle Counting  Models");

    println!("Crack growth is calculated for each cycle. The cyles can
be input directly or extracted from a sequence. The way the cycles are
determined affects the growth. The methods for extracting cycles from
a sequence are:");

    println!(
        "
rainflow  Crack growth is calculated from rainflow cycles i.e. the
          stress intensity range comes from the range of the rainflow
          cycles. Note this has a slight re-ordering effect that may upset
          the order of any image plots created.

tension   Crack growth calculated from tension part of cycle i.e. from a valley
          to the next peak. The striation pattern follows these tension cycles.
");

    header.section("da/dN data");
    println!(
        "The da/dN model consists of EQUATION:material where the
equation variable specifies the name of the da/dN equation and is
one of  [nasgro, paris, forman, walker, burchill, hartman, white, file]
The material variable specifies the name of the parameters to use for
that equation. If the values are given in --parameters they will
be used instead of the standard library values.\n
"
    );

    println!("{:35} {:20} Coefficients", "Name", "Ref.");
    let materials = material::get_all_dadns();
    for mat in materials.iter() {
        print!("{:35} {:20} ", mat.name, mat.cite);
        for num in &mat.eqn.variables() {
            print!("{:.5e} ", num)
        }
        println!();
    }
    println!(
        "{:35} {:20} {:30}",
        "file:FILE", " ", "Read FILE of tabular dadn data."
    );

    header.section("File formats");

    for &(file, form) in &formats {
        header.subsection(file);
        println!("{}", form);
    }

    header.section("References");

    for bib in &biblio {
        println!("{} {}\n", bib[0], bib[1]);
    }

    println!();
}

struct Counter {
    section: usize,
    subsection: usize,
}

impl Counter {
    fn new() -> Counter {
        Counter {
            section: 0,
            subsection: 0,
        }
    }

    // print as a header
    fn section(&mut self, head: &str) {
        self.section += 1;
        let header = format!("{}. {}", self.section, head);
        println!("\n{}", header);
        // Underline
        for _ in 0..header.len() {
            print!("=");
        }
        println!("\n");
    }

    fn subsection(&mut self, head: &str) {
        self.subsection += 1;
        let header = format!("{}.{}. {}", self.section, self.subsection, head);
        println!("{}", header);
        // Underline
        for _ in 0..header.len() {
            print!("-");
        }
        println!("\n");
    }
}
