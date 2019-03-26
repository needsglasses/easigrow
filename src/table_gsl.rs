extern crate rgsl; // used for spline interpolation of table
use log::debug;

use table;

/// increase the table so there is a row column for each column.
pub struct PairTable {
    /// columns
    pub columns: Vec<f64>,
    /// rows
    pub rows: Vec<Vec<f64>>,
    /// vectors of columns
    pub values: Vec<Vec<f64>>,
    // GSL Splines
    splines: Vec<rgsl::Spline>,
}

impl PairTable {
    pub fn new(columns: Vec<f64>, rows: Vec<Vec<f64>>, values: Vec<Vec<f64>>) -> PairTable {
        let mut splines = Vec::new();

        for i in 0..columns.len() {
            let spline = rgsl::Spline::new(&rgsl::InterpType::cspline(), rows[i].len()).unwrap();
            spline.init(&rows[i], &values[i]);
            splines.push(spline);
        }

        PairTable {
            columns,
            rows,
            values,
            splines,
        }
    }

    /// This is a simple interpolation of a table of variables.
    ///
    /// The most signicant variable (the one that produces the largest
    /// change in the return variable) should be the rows (which are
    /// spline fitted for accuracy) whereas the slower secondary
    /// variable are the coluumns (which are only linearly
    /// interpolate). Splines are fit to the upper and lower columns
    /// around the desired `column` point and then interpolate to the
    /// row value for each column. The lower and upper interpolations
    /// values are linearly interpolated to the column value.
    ///
    /// If `invert` is true we swap the `values` and `rows` so that we
    /// are interpolating the `values` to find the `row` like
    /// value. This routine does not extrapolate so both the row and
    /// column variables are capped to be the maximum or minimum values
    /// in the table.
    pub fn interp(&self, row: f64, column: f64) -> f64 {
        debug!("New interpolation row: {}, column: {}", row, column);
        let col_limited = column
            .max(self.columns[0])
            .min(self.columns[(self.columns.len() - 1)]);
        let col_low = table::nearest(col_limited, &self.columns);
        let row_limited = row.max(self.rows[col_low][0])
            .min(self.rows[col_low][(self.rows[col_low].len() - 1)]);
        debug!("Col limited {:?}, row limited {:?}", col_limited, row_limited );
        debug!("col {}", col_low);
        let mut accel = rgsl::InterpAccel::new();

        debug!("Columns: {:?}", self.columns);
        debug!("interp for {} using col {}", row, col_low);
        let row_s = self.splines[col_low].eval(row_limited, &mut accel);
        debug!("done interp row_s {}", row_s);

        // simple linear interpolation between the lower and upper columns
        if self.columns.len() > 1 {
            debug!("End interpolation");
            let mut accel = rgsl::InterpAccel::new();
            let row_limited = row.max(self.rows[col_low + 1][0])
                .min(self.rows[col_low + 1][(self.rows[col_low + 1].len() - 1)]);
            let row_e = self.splines[col_low + 1].eval(row_limited, &mut accel);
            debug!("done interp row_e {}", row_e);
            let interp = ((column - self.columns[col_low])
                / (self.columns[col_low + 1] - self.columns[col_low]))
                * (row_e - row_s) + row_s;
            debug!("interp result: {}", interp);
            interp
        } else {
            row_s
        }
    }
}
