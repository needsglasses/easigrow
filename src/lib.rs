//! fatigue
//! A library of useful fatigue functions
//! Paul White (2015)

pub static COMMENT: &'static str = "#  ";

extern crate svg;

#[cfg(not(feature = "GSL"))]
extern crate bspline;

#[cfg(not(feature = "GSL"))]
pub mod table_bspline;

#[cfg(feature = "GSL")]
pub mod table_gsl;

pub mod beta;
pub mod cycle;
pub mod dadn;
pub mod grow;
pub mod tag;
pub mod io;
pub mod plastic;
pub mod material;
pub mod fracto;
pub mod numbers;
pub mod table;
