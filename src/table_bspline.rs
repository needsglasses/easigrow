/// Tables with spline interpolation using the bspline library

use log::debug;

use bspline;
use table::nearest;

pub struct PairTable {
    /// columns
    pub columns: Vec<f64>,
    /// rows
    pub rows: Vec<Vec<f64>>,
    /// vectors of columns
    pub values: Vec<Vec<f64>>,
    /// Splines
    splines: Vec<bspline::BSpline<f64>>,
}

impl PairTable {
    pub fn new(columns: Vec<f64>,
               rows: Vec<Vec<f64>>,
               values: Vec<Vec<f64>>) -> PairTable {
        let mut splines = Vec::new();
        let degree = 2;
        debug!("PairTable::new()");
        debug!("columns: {:?}", columns);
        debug!("rows: {:?}", rows);
        debug!("values: {:?}", values);

        // make a spline for every set of row and column
        for i in 0..columns.len() {
            let mut knots = Vec::new();
            let start = rows[i][0];
            knots.push(start);
            knots.push(start);
            let mut t = rows[i].clone();
            knots.append(&mut t);
            let end = knots[knots.len() - 1];
            knots.push(end.clone());

            let mut points = Vec::new();
            let mut v = values[i].clone();

            points.append(&mut v);
            
            debug!("Knots are: {:?}", knots);
            let spline = bspline::BSpline::new(degree, points, knots);
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
        let col_low = nearest(col_limited, &self.columns);
        let row_limited = row.max(self.rows[col_low][0])
            .min(self.rows[col_low][(self.rows[col_low].len() - 1)]);
               debug!("Col limited {:?}, row limited {:?}", col_limited, row_limited );
               debug!("col {}", col_low);

               debug!("Columns: {:?}", self.columns);
               debug!("interp for {} using col {}", row, col_low);
               debug!("knot domain {:?}", self.splines[col_low].knot_domain());

        let row_s = self.splines[col_low].point(row_limited);
        debug!("done interp row_s {}", row_s);      
        debug!("No. of columns {}", self.columns.len());

        // simple linear interpolation between the lower and upper columns
        if self.columns.len() > 1 {
            debug!("End interpolation");
            let row_limited = row.max(self.rows[col_low + 1][0])
                .min(self.rows[col_low + 1][(self.rows[col_low + 1].len() - 1)]);
            let row_e = self.splines[col_low + 1].point(row_limited);
            debug!("done interp row_e {:?}", row_e);
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
