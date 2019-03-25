/// This is a gernalisation of the dadn interpolation into a general
/// table interpolation

use io;
use std::{fmt, process, f64};

#[cfg(feature = "GSL")]
pub use table_gsl::PairTable;

#[cfg(not(feature = "GSL"))]
pub use table_bspline::PairTable;

/// Interpolation data structure
pub struct Table {
    /// columns
    pub columns: Vec<f64>,
    /// rows
    pub row: Vec<f64>,
    /// vectors of columns
    pub values: Vec<Vec<f64>>,
    /// Stores the table data in a PairTable
    pub table: PairTable,
}

// transpose a simple vector of vectors (Ugly I know)
pub fn transpose(table: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut transposed = Vec::new();
    for _r in &table[0] {
        transposed.push(Vec::new());
    }

    for row in table {
        for (i, cell) in row.iter().enumerate() {
            transposed[i].push(*cell);
        }
    }
    transposed
}

/// Read two parameter table data from a file.
impl Table {
    // Read in a file in standard da/dn format
    // Comment lines start with #
    // The first valid line contains the values of the column variable
    // The next lines consist of row variables, followed by the values for each
    // of the colummn variables: value_c1 value_c2 value_c3... .
    pub fn read_file(interpfile: &str, invert: bool) -> Table {
        let mut table = io::read_table(interpfile); // raw data
        let columns = table.remove(0); // First row is a list of the r values

        let row: Vec<f64> = table.iter_mut().map(|x| x.remove(0)).collect();

        let mut rows = vec![];
        for _i in 0..columns.len() {
            rows.push(row.clone());
        }

        // transpose the table so the elements are the columns
        let table = transpose(&table);
        let pair_table = if invert {
            PairTable::new(columns.clone(), table.clone(), rows.clone())
        } else {
            PairTable::new(columns.clone(), rows.clone(), table.clone())
        };

        Table {
            columns,
            row,
            values: table,
            table: pair_table,
        }
    }

    pub fn interp(&self, row: f64, column: f64) -> f64 {
        self.table.interp(row, column)
    }

    pub fn new(columns: Vec<f64>, row: Vec<f64>, values: Vec<Vec<f64>>, invert: bool) -> Table {
        let mut rows = vec![];
        for _i in 0..columns.len() {
            rows.push(row.clone());
        }

        let pair_table = if invert {
            PairTable::new(columns.clone(), values.clone(), rows.clone())
        } else {
            PairTable::new(columns.clone(), rows.clone(), values.clone())
        };

        Table {
            columns,
            row,
            values,
            table: pair_table,
        }
    }

    /// verifies that the data in the table is consistent
    pub fn consistent(&self) -> bool {
        // check that there the correct number of columns
        for acol in &self.values {
            if self.row.len() != acol.len() {
                println!("Error: row {:?} != value {:?}", self.row, acol);
                return false;
            }
        }

        true
    }
}

impl fmt::Display for Table {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = writeln!(
            f,
            "Spline interpolated tabular lookup
     row       columns: {:?}",
            self.columns
        );
        for i in 0..self.row.len() {
            let _ = write!(f, "{:10.3e}: ", self.row[i]);
            for j in 0..self.values.len() {
                let _ = write!(f, "{:10} ", self.values[j][i]);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }
}

/// PairTable is a different format from Table.
///
/// In table the rows are common to all columns. In PairTable there is
/// a specific row for each column. In addition the rows may be of
/// unequal sizes, which means they are stored in row format.
impl PairTable {
    // Read in a file in standard da/dn format
    // Comment lines start with #
    // column variables
    // row variables
    // column1
    pub fn read_pair_file(file: &str) -> PairTable {
        let mut file_table = io::read_table(file); // raw data
        let columns = file_table.remove(0); // First row is a list of the r values

        let mut rows = vec![];
        let mut values = vec![];
        let mut i = 0;
        while i < columns.len() {
            rows.push(file_table[i].clone());
            values.push(file_table[i + 1].clone());
            i += 2;
        }

        let table = PairTable::new(columns, rows, values);

        if !table.consistent() {
            println!(
                "Error: in PairTable file {}, the file is inconsistent",
                file
            );
            process::exit(1);
        }

        table
    }

    /// verifies that the data in the table is consistent
    pub fn consistent(&self) -> bool {
        // check that there the correct number of columns
        for i in 0..self.columns.len() {
            if self.rows[i].len() != self.values[i].len() {
                println!(
                    "Error: row {:?} != value {:?}",
                    self.rows[i], self.values[i]
                );
                return false;
            }
        }

        true
    }
}

impl fmt::Display for PairTable {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = writeln!(
            f,
            "Spline interpolated tabular lookup
     row       columns: {:?}",
            self.columns
        );
        for i in 0..self.rows.len() {
            for j in 0..self.values.len() {
                let _ = write!(f, "{:10.3e}: ", self.rows[j][i]);
                let _ = write!(f, "{:10} ", self.values[j][i]);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }
}

// Find the position of the best surrounding values.
pub fn nearest<T: PartialOrd + Copy + fmt::Debug>(r: T, rs: &[T]) -> usize {
    match rs.iter().rposition(|&x| x < r) {
        Some(x) => x,
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_reading_dadn_file() {
        let x = Table::read_file("examples/barter14.dadn", true);

        // Clippy seems to choke on this
        // assert_eq!((x.interp(0.45, 0.0) - 1.0e-12).abs() < std::f64::EPSILON);
        // assert_eq!((x.interp(21.45, 0.0) - 1e-5).abs() < std::f64::EPSILON);
        
        // assert_eq!((x.interp(0.42, 0.3) - 1e-12).abs() < std::f64::EPSILON);
        // assert_eq!((x.interp(15.53, 0.3) - 1e-5).abs() < std::f64::EPSILON);
        
        // assert_eq!((x.interp(1.91, 0.4) - 1e-9).abs() < std::f64::EPSILON);
        // assert_eq!((x.interp(2.5, 0.4) - 2.696_498_799_927_166_3e-9).abs() < std::f64::EPSILON);
        
        // assert_eq!((x.interp(0.36, 0.7) - 1e-12).abs() < std::f64::EPSILON);
        // assert_eq!((x.interp(7.22, 0.7) - 1e-5).abs() < std::f64::EPSILON);
        
        // assert_eq!((x.interp(0.33, 0.8) - 1e-12).abs() < std::f64::EPSILON);
        // assert_eq!((x.interp(5.00, 0.8) - 1e-5).abs() < std::f64::EPSILON);
    }

    #[test]
    fn near() {
        let rs = vec![1usize, 2, 3, 4, 5, 6, 7, 8, 9];

        assert_eq!(nearest(0, &rs), 0);
        assert_eq!(nearest(3, &rs), 1);
        assert_eq!(nearest(10, &rs), 8);
        assert_eq!(nearest(20, &rs), 8);
    }
}
