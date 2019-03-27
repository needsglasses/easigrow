//! Routines for reading and writing files.

use std::io::{Read, Write};
use std::fs::File;
use tag::Tag;
use cycle::Cycle;
use std::path::Path;
use std;
use log::error;

use std::error::Error;

/// structure to hold measurements
#[derive(Debug, Clone)]
pub struct Measurement {
    pub line: usize,
    pub block: f64,
    pub a: f64,
}

/// Read in crack growth measurements from a file.
///
/// The crack data is given in a file containing:
/// '<lineno> <block> a'
///
/// a blank line indicates a break in the continuity of the
/// measurements. <lineno> is the corresponding line in the load file
/// that the measurement is made (starting at 0). If the line number
/// is not given it is assumed to be 0.  Unless there are measurements
/// made at smaller than 1 block intervals the assumption of 0 will
/// not affect the predictions.  If <block> is not given then it is
/// assumed the measurements are 1 block apart.
///
/// For example:
///  0 2 0.010
/// 15 2 0.013
///  0 3 0.015
///
///  0 5 0.020
/// 15 5 0.024
///  0 6 0.028
pub fn read_fracto_file(crack_file: &str, nseq: usize) -> Vec<Measurement> {
    let meas = read_table(crack_file);
    let mut fracto = Vec::new();

    for (b, m) in meas.iter().enumerate() {
        let f = match m.len() {
            0 => Measurement {
                line: 0,
                block: 0.0,
                a: 0.0,
            },
            1 => Measurement {
                line: 0,
                block: b as f64,
                a: m[0],
            },
            2 => Measurement {
                line: 0,
                block: m[0],
                a: m[1],
            },
            3 => Measurement {
                line: m[0] as usize,
                block: m[1] + m[0] / nseq as f64,
                a: m[2],
            },
            _ => {
                error!("Error: unknown format for fracto file '{}'.", crack_file);
                std::process::exit(1)
            }
        };
        fracto.push(f);
    }

    fracto
}

/// Read a table of data from a file.
pub fn read_table(filename: &str) -> Vec<Vec<f64>> {
    let path = Path::new(filename);
    let display = path.display();

    let mut file = match File::open(&path) {
        // The `description` method of `io::Error` returns a string that
        // describes the error
        Err(why) => {
            error!(
                "Error: could not open the file '{}': {}.",
                display,
                Error::description(&why)
            );
            std::process::exit(1)
        }
        Ok(file) => file,
    };

    let mut s = String::new();
    match file.read_to_string(&mut s) {
        Err(why) => {
            error!(
                "Error: could not read {}: {}.",
                display,
                Error::description(&why)
            );
            std::process::exit(1)
        }
        Ok(_) => true,
    };

    let mut arr = Vec::with_capacity(s.len() / 5); // make a guess at the capacity by assuming 5 characters per line
    for l in s.lines() {
        // skip any lines starting with comment character
        match l.chars().next() {
            Some(char) => if char == '#' {
                continue;
            },
            None => (),
        }

        let nums: Vec<f64> = l.split_whitespace()
            .map(|number| number.parse().unwrap())
            .collect();
        arr.push(nums);
    }

    arr
}

/// Read in a sequence.
///
/// There can be any number of values on a line in the file
/// and they will all be combined to return the sequence.
pub fn read_sequence(filename: &str) -> Vec<Tag> {
    let file_seq = read_table(filename);
    let mut seq = Vec::with_capacity(file_seq.len());

    for v in file_seq {
        for s in v {
            seq.push(s);
        }
    }
    Tag::from(&seq)
}

/// Write out a simple sequence to a file.
pub fn write_sequence(filename: &str, seq: &[Tag]) {
    let path = Path::new(filename);

    let mut file = match File::create(&path) {
        Err(why) => {
            error!("Error: Could not create file '{}': {}.", filename, &why);
            std::process::exit(1)
        }
        Ok(file) => file,
    };

    for load in seq.iter() {
        let result = writeln!(&mut file, "{}", load.value);
        match result {
            Ok(_result) => (),
            Err(err) => error!(
                "Error: unable to sucessfully write the sequence file - '{}'",
                err
            ),
        }
    }
}

/// Read in a cycle file in AFGROW format: 'max min n'
///
/// Where the number of cycles is greater than 1 we simply repeat the
/// cycles n in the cycle list to allow us to grow the sequence a
/// cycle at a time.
pub fn read_afgrow_cycles(cyclesfile: &str) -> Vec<Cycle<Tag>> {
    let mut rcycles = read_table(cyclesfile);
    let mut cycles = Vec::new();
    let mut count = 0;
    rcycles.remove(0); // skip the first row since it is the 'nsubspectrum nlevels' required by AFGROW

    for c in rcycles {
        for _ in 0..c[2] as usize {
            cycles.push(Cycle {
                max: Tag::new(c[0], count),
                min: Tag::new(c[1], count + 1),
            });
            count += 2;
        }
    }

    cycles
}

/// Write out the cycles in AFGROW format: 'max min 1'
/// i.e. each cycle is written out on a separate line.
pub fn write_cycles(file: &str, cycles: &[Cycle<Tag>]) {
    let mut file = match File::create(&file) {
        Err(why) => {
            error!("Error: Could not create file '{:?}': {}.", file, why);
            std::process::exit(1)
        }
        Ok(file) => file,
    };

    writeln!(file, "0 {}\n", cycles.len()).unwrap();

    for &Cycle { min, max } in cycles.iter() {
        let result = writeln!(&mut file, "{:?} {:?} 1\n", max.value, min.value);
        match result {
            Ok(_result) => (),
            Err(err) => error!(
                "Error: unable to sucessfully write the cycle file '{}'",
                err
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate tempdir;

    use self::tempdir::TempDir;
    use cycle;
    use tag::Tag;
    use io::*;

    #[test]
    fn write_cycle_file() {
        let tmp_dir = TempDir::new("testing").unwrap();
        let path = tmp_dir.path().join("test-cycle.out");

        let seq = vec![0.0, 3.0, 2.0, 5.0, 0.0];
        let (rain, _) = cycle::rainflow(&Tag::from(&seq));

        write_cycles(&path.to_str().unwrap(), &rain);
    }
}
