//! A variety of beta solutions for cracks in different shaped coupons.
//!
//! Beta functions relate the far field loading to the crack tip so
//! they are a function of the geometry and of the size of the crack
//! in the component. In all the functions here we use
//! non-dimensionalised factors. All the betas can be called with the
//! same set of non-dimensionalised parameters about the crack size
//! but each particular functino will only use those vaiables that are
//! relavent to it.
//!
//! The naming convention is (crack shape, geometry, loading) - (author, publication year)
//! The letter for crack is not included in the abbreviation since it seems redundant.
//! No need to say surface for semi-elliptical or quarter cracks.

// cargo test -- --nocapture

use table;
use std::process;
use std::f64::consts::{FRAC_PI_2, PI};
use grow;
use std::fmt;

use log::{info, error};

// This trait is constructed for betas to allow extra information to
// be supplied with the betas such as required for a tabular beta or
// the coupon beta which uses the dimensions of the coupon.
pub trait Beta {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, a_on_r: f64, phis: &[f64]) -> Vec<f64>;
    fn area(&self, _a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0
    }

    /// Convert a beta function into a beta table
    fn as_table(&self) -> TableBeta {
        let a_on_cs = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let a_on_ds = vec![
            0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13,
            0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25,
        ];

        let phis = vec![0.0, FRAC_PI_2];

        let mut values = Vec::new();
        for a_on_c in &a_on_cs {
            let mut column = Vec::new();
            for a_on_d in &a_on_ds {
                column.push(self.beta(*a_on_d, *a_on_c, 0.0, 0.0, &phis)[0]);
            }
            values.push(column);
        }

        let table = table::Table::new(a_on_cs, a_on_ds, values, false);
        TableBeta { table }
    }
}

impl fmt::Display for TableBeta {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = writeln!(f, "# {:9}", "a/d \\ a/c");
        let _ = write!(f, "            ");
        for a_on_c in &self.table.columns {
            let _ = write!(f, "{:6.3} ", a_on_c);
        }
        let _ = writeln!(f);

        for i in 0..self.table.row.len() {
            let _ = write!(f, "{:9.3}   ", self.table.row[i]);
            for j in 0..self.table.columns.len() {
                let value = self.table.values[j][i];
                let _ = write!(f, "{:6.3} ", value);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }
}

/// This is a generic beta value
/// Some of these betas do not use the phi vector and so the consistentancy is wrong.
pub fn get_beta_fn(beta_name: &str, component: &grow::Component) -> Box<Beta> {
    let all_betas = get_all_betas(component);

    if beta_name.starts_with("file:") {
        let components: Vec<&str> = beta_name.split(':').collect();
        let table = table::Table::read_file(components[1], false);
        Box::new(TableBeta { table })
    } else {
        match all_betas.into_iter().find(|x| x.name == beta_name) {
            Some(beta) => beta.eqn,
            None => {
                error!("Error: Unknown beta type {}", beta_name);
                process::exit(1);
            }
        }
    }
}

pub struct BetaCite<'a> {
    pub name: &'a str,
    pub summary: &'a str,
    pub cite: &'a str,
    pub args: &'a str,
    pub eqn: Box<Beta>,
}

// all semi-elliptical and quarter  cracks are surface cracks
pub fn get_all_betas(component: &grow::Component) -> Vec<BetaCite<'static>> {
    vec![
        BetaCite {
            name: "qct-broek86",
            summary: "quarter circular crack in an infinite plate in tension",
            cite: "[broek86]",
            args: "",
            eqn: Box::new(QuarterBroek86 {}),
        },
        BetaCite {
            name: "seft-newman84",
            summary: "semi-elliptical surface crack in a finite plate in tension",
            cite: "[Newman79]",
            args: "a/d, a/c, c/b, phi",
            eqn: Box::new(SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {}),
        },
        BetaCite {
            name: "seit-anderson05",
            summary: "semi-elliptical surface crack in an infinite plate in tension",
            cite: "[Anderson05]",
            args: "a/c, phi",
            eqn: Box::new(SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {}),
        },
        BetaCite {
            name: "qcft-murakami87",
            summary: "quarter circular corner crack in a finite plate in tension",
            cite: "[Murakami87]",
            args: "a/d",
            eqn: Box::new(QuarterCircularCornerCrackFinitePlateTensionMurakami87 {}),
        },
        BetaCite {
            name: "qeft-newman84",
            summary: "quarter elliptical corner crack in a finite plate in tension",
            cite: "[Newman79]",
            args: "a/d, a/c, c/b, phi",
            eqn: Box::new(QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {}),
        },
        BetaCite {
            name: "eft-newman84",
            summary: "elliptical crack in a finite plate in tension",
            cite: "[Newman79]",
            args: "a/d, a/c, c/b, phi",
            eqn: Box::new(EllipticalEmbeddedCrackFinitePlateTensionNewman84 {}),
        },
        BetaCite {
            name: "sset-tada73",
            summary: "single sided edge crack in a plate in tension",
            cite: "[Tada73]",
            args: "a/d",
            eqn: Box::new(SingleSidedEdgeCrackTensionTada73 {}),
        },
        BetaCite {
            name: "dset-tada73",
            summary: "double sided edge crack in a plate in tension",
            cite: "[Tada73]",
            args: "a/d",
            eqn: Box::new(DoubleSidedEdgeCrackTensionTada73 {}),
        },
        BetaCite {
            name: "compact-tada73",
            summary: "compact specimen in tension (scale is in load units not stress units) ",
            cite: "[Tada73]",
            args: "a/d, depth, width",
            eqn: Box::new(CompactCoupon::new(component.sideways, component.forward)),
        },
        BetaCite {
            name: "ct-fedderson66",
            summary: "centre cracked plate in tension",
            cite: "[Fedderson66]",
            args: "a/d",
            eqn: Box::new(CentreCrackTensionFedderson66 {}),
        },
        BetaCite {
            name: "ct-koiter65",
            summary: "centre cracked plate in tension",
            cite: "[Koiter65]",
            args: "a/d",
            eqn: Box::new(CentreCrackTensionKoiter65 {}),
        },
        BetaCite {
            name: "qcct-mcdonald07",
            summary: "vertically constrained coupon with corner crack in tension",
            cite: "[McDonald07]",
            args: "a/d",
            eqn: Box::new(CornerCrackConstrainedTensionMcdonald07::new()),
        },
        BetaCite {
            name: "ssht-bowie56",
            summary: "single sided through crack in a circular hole in tension",
            cite: "[Bowie56]",
            args: "a/r",
            eqn: Box::new(SingleSidedThroughCrackCircularHoleTensionBowie56 {}),
        },
        BetaCite {
            name: "dsht-bowie56",
            summary: "double sided crack through in a circular hole in tension",
            cite: "[Bowie56]",
            args: "a/r",
            eqn: Box::new(DoubleSidedThroughCrackCircularHoleTensionBowie56 {}),
        },
        BetaCite {
            name: "dccht-newman81",
            summary: "double sided corner crack in a hole in tension",
            cite: "[Newman81]",
            args: "a/d, a/c, c/b, a/r, phi",
            eqn: Box::new(DoubleSidedCornerCrackHoleTensionNewman81 {}),
        },
        BetaCite {
            name: "serbb-shin04",
            summary: "semi-elliptical surface crack in a round bar in bending",
            cite: "[shin04]",
            args: "a/d, a/c",
            eqn: Box::new(SemiEllipticalSurfaceCrackRoundBarBendingShin04::new()),
        },
        BetaCite {
            name: "serbb-murakami87",
            summary: "semi-elliptical surface crack in a round bar in bending",
            cite: "[Murakami87]",
            args: "a/d, a/c",
            eqn: Box::new(SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {}),
        },
        BetaCite {
            name: "serbt-murakami87",
            summary: "semi-elliptical surface crack in a round bar in tension",
            cite: "[Murakami87]",
            args: "a/d, a/c",
            eqn: Box::new(SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {}),
        },
        BetaCite {
            name: "serbb-murakami86",
            summary: "semi-elliptical surface crack in a round bar in bending",
            cite: "[Murakami86]",
            args: "a/d, a/c",
            eqn: Box::new(SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {}),
        },
        BetaCite {
            name: "serbt-murakami86",
            summary: "semi-elliptical surface crack in a round bar in tension",
            cite: "[Murakami86]",
            args: "a/d, a/c",
            eqn: Box::new(SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {}),
        },
        BetaCite {
            name: "esb-murakami87",
            summary: "edge crack in a strip in bending",
            cite: "[Murakami87]",
            args: "a/d",
            eqn: Box::new(EdgeCrackStripBendingMurakami87 {}),
        },
        BetaCite {
            name: "est-murakami87",
            summary: "edge crack in a strip in tension",
            cite: "[Murakami87]",
            args: "a/d",
            eqn: Box::new(EdgeCrackStripTensionMurakami87 {}),
        },
        BetaCite {
            name: "file:FILE",
            summary: "read FILE for beta values",
            cite: "",
            args: "a/d, a/c",
            // dummy placeholder
            eqn: Box::new(EdgeCrackStripTensionMurakami87 {}),
        },
    ]
}

pub struct TableBeta {
    table: table::Table,
}

impl Beta for TableBeta {
    /// Convert a beta function into a beta table
    fn as_table(&self) -> TableBeta {
        let table = table::Table::new(
            self.table.columns.clone(),
            self.table.row.clone(),
            self.table.values.clone(),
            false,
        );
        TableBeta { table }
    }

    fn beta(
        &self,
        a_on_d: f64,
        a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        vec![self.table.interp(a_on_d, a_on_c)]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - a_on_d
    }
}

struct QuarterBroek86 {}

impl Beta for QuarterBroek86 {
    fn beta(
        &self,
        _a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        let embedded_beta = 2.0 / PI;
        vec![1.12 * 1.12 * embedded_beta]
    }
}

/// Beta factor for a corner crack in a plate in tension.
/// Plate is of thickness theta = 45 degrees.
///
/// Ref. Stress intensity factors handbook. Vol 2,
/// by Y. Murakami,
/// Pergamon Press, Oxford, , 1987.
struct QuarterCircularCornerCrackFinitePlateTensionMurakami87 {}

impl Beta for QuarterCircularCornerCrackFinitePlateTensionMurakami87 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        let a_on_t = a_on_d;
        vec![
            ((1.05 + (-0.44 + 1.06 / 1.3) * a_on_t.powi(2) - 0.25 * a_on_t.powi(4))
                * (1.0 + (0.08 + 0.40 * a_on_t.powi(2)) * (1.0 - (PI / 4.0).sin()).powi(3))
                * (1.0 + (0.08 + 0.15 * a_on_t.powi(2)) * (1.0 - (PI / 4.0).cos()).powi(3))
                * ((PI / 4.0).cos().powi(2) + (PI / 4.0).sin().powi(2)).powf(0.25)) / 2.464f64.sqrt()
        ]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        PI * a_on_d.powi(2) / 4.0
    }
}

/// Semi-elliptical beta factor for an edge crack.
/// Where *a* is the depth of crack, *2c*  is the surface length of the crack,
/// *`a_on_c` = a/c* and `phi` are the angle to a point on the crack front from surface .
///
/// Ref. Fracture Mechanics p.48,
/// by T. L. Anderson,
/// Taylor and Francis Group 2005.
///
/// Calculates the beta factor for a semi-elliptical crack in infinite
/// plate phi is the angle of the position on the crack front.  This
/// has been re-defined phi from the book so that it is from the
/// centreline i.e. phi = 0 is now at the deepest point and PI/2 is at
/// the surface.
struct SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {}

impl Beta for SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {
    fn beta(
        &self,
        _a_on_d: f64,
        a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        phis: &[f64],
    ) -> Vec<f64> {
        let mut beta = Vec::new();

        if a_on_c > 1.0 {
            println!(
                "Error: semi_elliptical_infinite_anderson05: invalid c > a. a/c = {}",
                a_on_c
            );
            process::exit(1);
        }

        for phi in phis {
            // surface correction factor
            let lambda_s = (1.13 - 0.09 * a_on_c) * (1.0 + 0.1 * (1.0 - phi.cos()).powi(2));

            // angle correction factor
            let f_phi = (phi.cos().powi(2) + (a_on_c * phi.sin()).powi(2)).powf(0.25);
            let q = 1.0 + 1.464 * a_on_c.powf(1.65);

            beta.push(lambda_s * f_phi / q.sqrt());
        }
        beta
    }
}

/// Beta factor for a semi-elliptical crack under tension.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// NASA Technical Memorandum 85793.
struct SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {}

impl Beta for SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, _a_on_r: f64, phis: &[f64]) -> Vec<f64> {
        let c_on_a = 1.0 / a_on_c;
        let mut beta = Vec::new();
        //    let mut h = 0.0f64;

        for phi in phis {
            let newman_phi = FRAC_PI_2 - phi;

            // finite width correction
            // sec x = 1/cos x
            let f_w = (FRAC_PI_2 * c_on_b * a_on_d.sqrt()).cos().recip().sqrt();

            let f = if a_on_c < 1.0 {
                let m1 = 1.13 - 0.09 * a_on_c;
                let m2 = -0.54 + 0.89 / (0.2 + a_on_c);
                let m3 = 0.5 - 1.0 / (0.65 + a_on_c) + 14.0 * (1.0 - a_on_c).powf(24.0);
                let g = 1.0 + (0.1 + 0.35 * a_on_d.powi(2)) * (1.0 - newman_phi.sin()).powi(2);
                info!("m1 {} m2 {} m3 {} g {} f_w {} f_phi {}", m1, m2, m3, g, f_w, f_phi(a_on_c, newman_phi));

                (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powi(4)) * g * f_w
                    * f_phi(a_on_c, newman_phi)
            } else {
                let m1 = c_on_a.sqrt() * (1.0 + 0.04 * c_on_a);
                let m2 = 0.2 * c_on_a.powf(4.0);
                let m3 = -0.11 * c_on_a.powf(4.0);

                // check on this it seems unsymmetrical
                // corrected the angle
                let g = 1.0
                    + (0.1 + 0.35 * c_on_a * a_on_d.powi(2))
                        * (1.0 - (FRAC_PI_2 - phi).sin()).powi(2);

                (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powi(4)) * g * f_w
                    * f_phi(a_on_c, newman_phi)
            };

            beta.push(f / shape_factor_newman84(a_on_c).sqrt())
        }
        beta
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        PI * a_on_d * c_on_b / 2.0
    }
}

/// Shape correction factor.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// Nasa Technical Memorandum 85793.

fn shape_factor_newman84(a_on_c: f64) -> f64 {
    if a_on_c < 1.0 {
        1.0 + 1.464 * a_on_c.powf(1.65)
    } else {
        1.0 + 1.464 * a_on_c.recip().powf(1.65)
    }
}

/// Beta factor for a quarter elliptical crack in tension.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// Nasa Technical Memorandum 85793.
struct QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {}

impl Beta for QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, _a_on_r: f64, phis: &[f64]) -> Vec<f64> {
        let c_on_a = 1.0 / a_on_c;
        let c_on_d = a_on_d / a_on_c;

        let mut beta = Vec::new();

        for phi in phis {
            let newman_phi = FRAC_PI_2 - phi;

            // finite width correction
            // let f_w = (PI * c_on_b * a_on_d.sqrt() /2.0).cos().sqrt().recip(); // murakami
            let lambda = c_on_b * a_on_d.sqrt();
            let f_w = 1.0 - 0.2 * lambda + 9.4 * lambda.powi(2) - 19.4 * lambda.powi(3)
                + 27.1 * lambda.powi(4);

            let f = if a_on_c <= 1.0 {
                let m1 = 1.08 - 0.03 * a_on_c;
                let m2 = -0.44 + 1.06 / (0.3 + a_on_c);
                let m3 = -0.5 + 0.25 * a_on_c + 14.8 * (1.0 - a_on_c).powf(15.0);
                let g1 = 1.0 + (0.08 + 0.4 * a_on_d.powi(2)) * (1.0 - newman_phi.sin()).powi(3);
                let g2 = 1.0 + (0.08 + 0.15 * a_on_d.powi(2)) * (1.0 - newman_phi.cos()).powi(3);

                // println!("a/c <= 1 m1 {} m2 {} m3 {} g1 {} g2 {}", m1, m2, m3, g1, g2 );
                (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g1 * g2 * f_w
                    * f_phi(a_on_c, newman_phi)
            } else {
                let m1 = c_on_a.sqrt() * (1.08 - 0.03 * c_on_a);
                let m2 = 0.375 * c_on_a.powi(2);
                let m3 = -0.25 * c_on_a.powi(2);

                // the report says c/t instead of c/b
                //            let g1 = 1.0 + (0.08 + 0.4 * c_on_b.powi(2)) * (1.0 - (FRAC_PI_2 - phi).sin()).powi(3);
                //            let g2 = 1.0 + (0.08 + 0.15 * c_on_b.powi(2)) * (1.0 - (FRAC_PI_2 - phi).cos()).powi(3);
                let g1 = 1.0 + (0.08 + 0.4 * c_on_d.powi(2)) * (1.0 - newman_phi.sin()).powi(3);
                let g2 = 1.0 + (0.08 + 0.15 * c_on_d.powi(2)) * (1.0 - newman_phi.cos()).powi(3);

                (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g1 * g2 * f_w
                    * f_phi(a_on_c, newman_phi)
            };

            // println!("a_on_c {} f_w {} f {} shape factor {} beta {}", a_on_c, f_w, f, shape_factor(a_on_c).sqrt(), f/shape_factor(a_on_c).sqrt());
            beta.push(f / shape_factor_newman84(a_on_c).sqrt())
        }
        beta
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        PI * a_on_d * c_on_b / 4.0
    }
}

/// Angular correction `f_phi`.
fn f_phi(a_on_c: f64, newman_phi: f64) -> f64 {
    if a_on_c <= 1.0 {
        (a_on_c.powi(2) * newman_phi.cos().powi(2) + newman_phi.sin().powi(2)).powf(0.25)
    } else {
        (a_on_c.recip().powi(2) * newman_phi.sin().powi(2) + newman_phi.cos().powi(2)).powf(0.25)
    }
}

/// Shape correction factor.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// NASA Technical Memorandum 85793.
struct EllipticalEmbeddedCrackFinitePlateTensionNewman84 {}

impl Beta for EllipticalEmbeddedCrackFinitePlateTensionNewman84 {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, _a_on_r: f64, phis: &[f64]) -> Vec<f64> {
        let c_on_a = 1.0 / a_on_c;
        let mut beta = Vec::new();

        for phi in phis {
            let newman_phi = FRAC_PI_2 - phi;
            // finite width correction
            let f_w = (FRAC_PI_2 * c_on_b * a_on_d.sqrt()).cos().sqrt();

            let m2 = 0.05 / (0.11 + a_on_c.powf(1.5));
            let m3 = 0.29 / (0.23 + a_on_c.powf(1.5));
            let g = 1.0
                - (a_on_d.powi(4) * (2.6 - 2.0 * a_on_d).sqrt() / (1.0 + 4.0 * a_on_c))
                    * newman_phi.cos().abs();

            let m1 = if a_on_c <= 1.0 { 1.0 } else { c_on_a.sqrt() };

            let f = (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g * f_w
                * f_phi(a_on_c, newman_phi);

            beta.push(f / shape_factor_newman84(a_on_c).sqrt())
        }

        beta
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        PI * a_on_d * c_on_b
    }
}

/// Double corner crack in a hole in tension.
///
/// Taken from Y. Murakami
/// Stress Intensity Factors handbook
/// Vol. 2. P. 716 Pergamon Press, 1987.
/// In turn taken from J. C. Newman Jr. and I. S. Raju
/// Stress Intensity Factor Equations for Cracks
/// in three-dimensional finite bodies, Nasa Technical Memorandum 83299 1981 p 1--49.
struct DoubleSidedCornerCrackHoleTensionNewman81 {}

impl Beta for DoubleSidedCornerCrackHoleTensionNewman81 {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, a_on_r: f64, phis: &[f64]) -> Vec<f64> {
        // println!("a_on_r {}, a_on_c {}", a_on_r, a_on_c);
        let c_on_a = 1.0 / a_on_c;
        let mut beta = Vec::new();
        let c_on_r = a_on_r / a_on_c;

        let m1 = if a_on_c <= 1.0 {
            1.13 - 0.09 * a_on_c
        } else {
            c_on_a.sqrt() * (1.0 + 0.04 * c_on_a)
        };
        let m2 = if a_on_c <= 1.0 {
            -0.54 + (0.89 / (0.2 + a_on_c))
        } else {
            0.2 * c_on_a.powf(1.5)
        };

        let m3 = if a_on_c <= 1.0 {
            0.5 - 1.0 / (0.65 + a_on_c)
        } else {
            -0.11 * c_on_a.powi(4)
        };

        let f_w = (FRAC_PI_2 * c_on_b * a_on_d.sqrt()).cos().sqrt();

        for phi in phis {
            let newman_phi = FRAC_PI_2 - phi;

            // finite width correction
            let g1 = if a_on_c < 1.0 {
                1.0 + (0.1 + 0.35 * c_on_a * a_on_d.powi(2)) * (1.0 - phi.sin()).powi(2)
            } else {
                1.0 + (0.1 + 0.35 * a_on_d.powi(2)) * (1.0 - phi.sin()).powi(2)
            };
            // println!("c_on_r {}", c_on_r);
            let lambda = 1.0 / (1.0 + c_on_r * (0.85 * phi).cos());
            let g2 = (1.0 - 0.15 * lambda + 3.46 * lambda.powi(2) - 4.47 * lambda.powi(3)
                + 3.52 * lambda.powi(4)) / (1.0 + 0.13 * lambda.powi(2));

            let g3 = if a_on_c <= 1.0 {
                (1.0 + 0.04 * a_on_c) * (1.0 + 0.1 * (1.0 - phi.cos()).powi(2))
                    * (0.8 + 0.2 * a_on_d.powf(0.25))
            } else {
                (1.13 - 0.09 * c_on_a) * (1.0 + 0.1 * (1.0 - phi.cos()).powi(2))
                    * (0.8 + 0.2 * a_on_d.powf(0.25))
            };

            let f = (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g1 * g2 * g3 * f_w
                * f_phi(a_on_c, newman_phi);
            // println!("m1 {}, m2 {}, m3 {}, lamda {}, g1 {}, g2 {}, g3 {}, f_w {}, f {}", m1, m2, m3, lambda, g1, g2, g3, f_w, f);
            beta.push(f)
        }

        beta
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, a_on_r: f64) -> f64 {
        1.0 - 2.0 * a_on_r - PI * a_on_d * c_on_b
    }
}

/// Single edge notched tension (SENT).
///
/// Ref. H. Tada, P.C. Paris and G. R. Irwin
/// The stress analysis of cracks handbook
/// P. 2.11
/// Compiled from Tada 1973
/// `a_on_d` is the ratio of crack depth a to specimen depth d
/// Uses the third equation which has accuracy better than 0.5% for any a/d
struct SingleSidedEdgeCrackTensionTada73 {}

impl Beta for SingleSidedEdgeCrackTensionTada73 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        let limited_a_on_d = a_on_d.min(0.8);
        let num = if limited_a_on_d == 0.0 {
            1.0
        } else {
            ((1.0 / (FRAC_PI_2 * limited_a_on_d)) * (FRAC_PI_2 * limited_a_on_d).tan()).sqrt()
        };

        let denom = (PI * limited_a_on_d / 2.0).cos();

        vec![
            (num / denom) * (0.752 + 2.02 * limited_a_on_d
                             + 0.37 * (1.0 - (FRAC_PI_2 * limited_a_on_d).sin()).powi(3))
        ]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - a_on_d
    }
}

/// Middle cracked tension (centre cracked tension 2nd eq)
///
/// Fedderson 1966
/// Taken from Damage Tolerant Design handbook from Afgrow documentation.
/// <http://www.afgrow.net/applications/DTDHandbook/Sections/page11_3.aspx#standard_center_cracked_tension_specimen>
/// Note that `a_on_d` divided by 2 to account for the different definition of the width/depth in the equation.
struct CentreCrackTensionFedderson66 {}

impl Beta for CentreCrackTensionFedderson66 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        let f_0 = (PI * a_on_d / 2.0).cos().powf(-0.5);

        vec![f_0]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - a_on_d
    }
}

/// Centre cracked tension.
///
/// Modification of Koiter 1965 from Tada 73
/// taken from H. Tada, P.C. Paris and G. R. Irwin (p 2.2)
struct CentreCrackTensionKoiter65 {}

impl Beta for CentreCrackTensionKoiter65 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        let denom = (1.0 - a_on_d).sqrt();
        let numerator = 1.0 - 0.5 * a_on_d + 0.370 * a_on_d.powi(2) - 0.044 * a_on_d.powi(3);

        vec![numerator / denom]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - a_on_d
    }
}

/// Double edge notched tension - DENT
///
/// Two edge cracks each of length a in a plate of width 2d
/// by H. Tada, P.C. Paris and G. R. Irwin.
struct DoubleSidedEdgeCrackTensionTada73 {}

impl Beta for DoubleSidedEdgeCrackTensionTada73 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        vec![
            (1.122 - 0.561 * a_on_d - 0.205 * a_on_d.powi(2) + 0.471 * a_on_d.powi(3)
                - 0.190 * a_on_d.powi(4)) / (1.0 - a_on_d).sqrt(),
        ]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - 2.0 * a_on_d
    }
}

/// Compact specimen width W, thickness B, 'a_on_d' is really a/w
///
/// by H. Tada, P.C. Paris and G. R. Irwin.
struct CompactCoupon {
    width: f64,
    depth: f64,
}

impl CompactCoupon {
    // compact_tension_tada73
    fn new(width: f64, depth: f64) -> CompactCoupon {
        CompactCoupon {
            width,
            depth,
        }
    }
}

impl Beta for CompactCoupon {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        // This beta factor is non-dimensionalised with load
        // hence we have to divide by  (PI * a_on_d * width).sqrt() so the will be cancelled out in the later calculation
        // this expression is valid for a/w > 0.2
        let correction = (PI * a_on_d * self.width).sqrt(); // = sqrt(PI * a)
        let front = (self.depth * self.width.sqrt()).recip();
        let scale = (2.0 + a_on_d) / (1.0 - a_on_d).powf(3.0 / 2.0);
        let rest = 0.886 + 4.64 * a_on_d - 13.32 * a_on_d.powi(2) + 14.72 * a_on_d.powi(3)
            - 5.60 * a_on_d.powi(4);

        let result = vec![front * scale * rest / correction];

        //    println!("a/d {} width {} thick {}", a_on_d, width, thick);
        //    println!("beta: a/d {} front {} scale {} rest {} correction {} result {:?} ", a_on_d, front, scale, rest, correction, result);
        result
    }
}

/// Semi-elliptical surface crack in a round bar in bending
/// Murakami87 p. 657
/// This is obtained by Murakami using the results for an edge crack in a strip
/// for which he has equations for both tension and bending whereas the crack in a
/// round bar is available for tension.
struct SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {}

impl Beta for SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, a_on_r: f64, phis: &[f64]) -> Vec<f64> {
        // println!("here in murakami bending");
        // println!("phis {:?}", phis);

        let semi_tension = SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {};
        let st = semi_tension.beta(a_on_d, a_on_c, c_on_b, a_on_r, phis);
        let edge_bending = EdgeCrackStripBendingMurakami87 {};
        let eb = edge_bending.beta(a_on_d, a_on_c, c_on_b, a_on_r, phis);
        let edge_tension = EdgeCrackStripTensionMurakami87 {};
        let et = edge_tension.beta(a_on_d, a_on_c, c_on_b, a_on_r, phis);

        // println!("st {:?}", st);
        // println!("eb {:?}", eb);
        // println!("et {:?}", et);
        
        let mut betas = Vec::new();
        for i in 0..phis.len() {
            betas.push(st[i] * eb[i] / et[i]);
        }

        betas
//        vec![st[0] * eb[0] / et[0], st[1] * eb[1] / et[1]]
    }
}

/// Stress Intensity Factor Equations for a Semi-Elliptical Surface Crack in a shaft under bending
/// Yukitaka Murakami and Hideto Tsuru 86
/// No. 87-0164B
/// This is obtained by Murakami using the results for an edge crack in a strip
/// for which he has equations for both tension and bending whereas the crack in a
/// round bar is available for tension.
struct SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {}

impl Beta for SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {
    fn beta(&self, a_on_d: f64, a_on_c: f64, c_on_b: f64, a_on_r: f64, phis: &[f64]) -> Vec<f64> {
        let semi_tension = SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {};
        let st = semi_tension.beta(a_on_d, a_on_c, c_on_b, a_on_r, phis);

        let edge_bending = EdgeCrackStripBendingMurakami87 {};
        let eb = edge_bending.beta(a_on_d, a_on_c, c_on_b, a_on_r, phis);

        let edge_tension = EdgeCrackStripTensionMurakami87 {};
        let et = edge_tension.beta(a_on_d, a_on_c, c_on_b, a_on_r, phis);

        let mut betas = Vec::new();
        for i in 0..phis.len() {
            betas.push(st[i] * eb[i] / et[i]);
        }

        betas
    }
}

/// Semi-elliptical surface crack in a round bar in tension
/// Murakami87 p. 657
struct SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {}

impl Beta for SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {
    // There seems to be an error somewhere as this equation does not
    // reproduce the table of sample points given by murakami87. It
    // comes close. The original reference provided is murakami and
    // Tsuru 86 which is also slightly different to murakami87 in the
    // second part of the equation and in the table of sample beta
    // factors for the round bar in bending.
    fn beta(
        &self,
        a_on_d: f64,
        a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {

        let beta = (1.122 - 0.230 * a_on_c - 0.901 * a_on_c.powi(2) + 0.949 * a_on_c.powi(3)
                    - 0.280 * a_on_c.powi(4))
            * (1.0 + 0.157 * a_on_d - 0.634 * a_on_d.powi(2) + 4.590 * a_on_d.powi(3)
                   - 6.628 * a_on_d.powi(4));

        vec![beta; _phis.len()]            
    }
}

/// Semi-elliptical surface crack in a round bar in tension
/// Murakami and Tsuru 86
struct SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {}

impl Beta for SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {
    fn beta(
        &self,
        a_on_d: f64,
        a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {

        // The original reference by
        // Stress Intensity Factor equations for a semi-elliptical surface crack in a shaft under bending
        let beta = (1.122 - 0.230 * a_on_c - 0.901 * a_on_c.powi(2) + 0.949 * a_on_c.powi(3) - 0.280 * a_on_c.powi(4)) *
            (1.0 + 0.314 * a_on_d - 2.536 * a_on_d.powi(2) + 36.72 * a_on_d.powi(3) - 106.048 * a_on_d.powi(4));

        vec![beta; _phis.len()]            
    }
}

/// Edge crack in a strup in bending
/// Murakami87 p. 657
struct EdgeCrackStripBendingMurakami87 {}

impl Beta for EdgeCrackStripBendingMurakami87 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {

        let beta = 1.121 - 1.199 * a_on_d + 4.775 * a_on_d.powi(2) - 1.628 * a_on_d.powi(3)
            - 7.035 * a_on_d.powi(4) + 13.27 * a_on_d.powi(5);
        
        vec![beta; _phis.len()]
    }
}

/// Edge crack in a strup in bending tension
/// Murakami87 p. 657
struct EdgeCrackStripTensionMurakami87 {}

impl Beta for EdgeCrackStripTensionMurakami87 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {

        let beta = 1.12 - 0.231 * a_on_d + 10.55 * a_on_d.powi(2) - 21.72 * a_on_d.powi(3)
            + 30.39 * a_on_d.powi(4);

        vec![beta; _phis.len()]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - a_on_d
    }
}

/// Beta for a circular corner/edge crack with the coupon constrained to extend uniaxially
///
/// The crack transitions from a corner to an edge crack
/// (Coupon is constrained so that ends are not free to rotate).
///
/// Here the a is in the depth direction of the coupon but it should
/// equal c The coupon is constrained so it is not free to rotate
/// which is more representative of a coupon clamped in the jaws in a
/// test machine. d is the thickness of the coupon.
/// This is probably for a 160 mm radius coupon but it is not clear.
///
/// Ref. Beta Values for Low Kt Specimens by M. `McDonald`
/// DSTO Minute Air07/048 Combat Aircraft Support.
struct CornerCrackConstrainedTensionMcdonald07 {
    table: table::Table,
}

impl CornerCrackConstrainedTensionMcdonald07 {
    fn new() -> TableBeta {
        let a = vec![
            0.0, 0.0001, 0.0006, 0.0011, 0.0016, 0.0021001, 0.0026003, 0.0031004, 0.0036005,
            0.0041007, 0.0046014, 0.005102, 0.0056027, 0.0061034, 0.0066051, 0.0071055, 0.0076062,
            0.0081089, 0.0086145, 0.0091147, 0.00962, 0.0101257,
        ];

        let betas = vec![
            0.709411963f64,
            0.709411963,
            0.710173817,
            0.714482118,
            0.727905309,
            0.750298769,
            0.778871362,
            0.808379238,
            0.845057485,
            0.892278507,
            0.952439289,
            1.013915829,
            1.076584454,
            1.134948417,
            1.188004837,
            1.240923778,
            1.293874445,
            1.347036619,
            1.374464709,
            1.405714819,
            1.434011108,
            1.46362585,
        ];

        // non-dimensionalise the crack depth data by the coupon width 25 mm
        let a_on_ds = a.iter().map(|x| x / 25.0e-3).collect::<Vec<f64>>();

        TableBeta {
            table: table::Table::new(vec![0.0], a_on_ds, vec![betas], false),
        }
    }
}

impl Beta for CornerCrackConstrainedTensionMcdonald07 {
    fn beta(
        &self,
        a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        vec![self.table.interp(a_on_d, 0.0)]
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0 - PI * a_on_d * a_on_d
    }
}

// Additional beta factors from
// @TechReport{paul88,
//   author = 	 {J. Paul and D. Lombardo},
//   title = 	 {CRKGRW - Crack growth program users manual},
//   institution =  {Defence Science and Technology Organisation},
//   year = 	 1988,
//   month = 	 {June}}

/// Bowie solution for a single crack from a circular hole
struct SingleSidedThroughCrackCircularHoleTensionBowie56 {}

impl Beta for SingleSidedThroughCrackCircularHoleTensionBowie56 {
    fn beta(
        &self,
        _a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        vec![0.6762 + 0.8734 / (0.3246 + a_on_r)]
    }
}

/// Bowie solution for a double cracked circular hole
struct DoubleSidedThroughCrackCircularHoleTensionBowie56 {}

impl Beta for DoubleSidedThroughCrackCircularHoleTensionBowie56 {
    fn beta(
        &self,
        _a_on_d: f64,
        _a_on_c: f64,
        _c_on_b: f64,
        a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        vec![0.9439 + 0.6865 / (0.2772 + a_on_r)]
    }
}

/// Another beta factor for a round bar with greater depth and hopefully more accurate than Murakami.
/// It provides the beta factor at two points around the crack front.
///
/// Experimental and finite element analyses on stress intensity
/// factors of an elliptical surface crack in a circular shaft under
/// tension and bending
/// C. S. Shin and C. Q. Cai
/// International Journal of Fracture 129: 239–264, 2004
struct SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
    deep: table::Table,
    surface: table::Table,
}

impl SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
    fn new() -> SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
        let a_on_ds = vec![
            0.067, 0.133, 0.200, 0.267, 0.333, 0.400, 0.467, 0.533, 0.600, 0.667, 0.733, 0.800
        ];
        let a_on_cs = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

        // x/h = 0.0
        let betas_deep = vec![
            vec![0.963, 0.954, 0.929, 0.878, 0.834, 0.786, 0.739, 0.692, 0.649, 0.609, 0.576],
            vec![0.897, 0.890, 0.870, 0.840, 0.801, 0.757, 0.710, 0.662, 0.618, 0.576, 0.537],
            vec![0.872, 0.866, 0.848, 0.820, 0.783, 0.739, 0.690, 0.640, 0.592, 0.547, 0.506],
            vec![0.879, 0.873, 0.856, 0.828, 0.790, 0.743, 0.692, 0.637, 0.583, 0.532, 0.486],
            vec![0.917, 0.911, 0.893, 0.863, 0.823, 0.773, 0.716, 0.654, 0.592, 0.532, 0.478],
            vec![0.991, 0.984, 0.964, 0.932, 0.888, 0.832, 0.767, 0.695, 0.621, 0.549, 0.482],
            vec![1.112, 1.104, 1.082, 1.045, 0.994, 0.930, 0.854, 0.768, 0.678, 0.588, 0.504],
            vec![1.302, 1.294, 1.268, 1.224, 1.164, 1.087, 0.995, 0.889, 0.775, 0.659, 0.550],
            vec![1.609, 1.599, 1.566, 1.512, 1.437, 1.341, 1.224, 1.088, 0.938, 0.783, 0.634],
            vec![2.126, 2.113, 2.070, 1.998, 1.899, 1.771, 1.614, 1.429, 1.222, 1.002, 0.787],
            vec![3.082, 3.063, 3.002, 2.899, 2.755, 2.570, 2.342, 2.069, 1.758, 1.421, 1.083],
            vec![5.140, 5.110, 5.011, 4.841, 4.606, 4.302, 3.923, 3.466, 2.934, 2.344, 1.737],
        ];

        // x/h = 1.0
        let betas_surface = vec![
            vec![0.486, 0.523, 0.553, 0.578, 0.596, 0.609, 0.616, 0.618, 0.613, 0.603, 0.587],
            vec![0.510, 0.548, 0.579, 0.604, 0.623, 0.635, 0.641, 0.640, 0.633, 0.619, 0.599],
            vec![0.557, 0.596, 0.629, 0.654, 0.673, 0.684, 0.689, 0.686, 0.677, 0.660, 0.637],
            vec![0.600, 0.640, 0.673, 0.699, 0.717, 0.728, 0.732, 0.728, 0.717, 0.699, 0.674],
            vec![0.654, 0.695, 0.729, 0.755, 0.774, 0.784, 0.788, 0.783, 0.771, 0.751, 0.724],
            vec![0.742, 0.786, 0.822, 0.850, 0.869, 0.880, 0.882, 0.877, 0.862, 0.840, 0.809],
            vec![0.877, 0.926, 0.966, 0.996, 1.017, 1.028, 1.029, 1.022, 1.004, 0.978, 0.941],
            vec![1.062, 1.118, 1.163, 1.197, 1.220, 1.231, 1.232, 1.222, 1.200, 1.168, 1.124],
            vec![1.326, 1.393, 1.446, 1.485, 1.511, 1.524, 1.524, 1.510, 1.483, 1.443, 1.389],
            vec![1.755, 1.838, 1.904, 1.953, 1.985, 2.000, 1.998, 1.979, 1.943, 1.890, 1.820],
            vec![2.544, 2.659, 2.749, 2.816, 2.859, 2.878, 2.873, 2.844, 2.791, 2.714, 2.614],
            vec![4.138, 4.317, 4.458, 4.560, 4.625, 4.651, 4.639, 4.589, 4.500, 4.374, 4.209],
        ];

        let deep = table::Table::new(
            a_on_cs.clone(),
            a_on_ds.clone(),
            table::transpose(&betas_deep),
            false,
        );
        let surface = table::Table::new(
            a_on_cs.clone(),
            a_on_ds.clone(),
            table::transpose(&betas_surface),
            false,
        );

        SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
            deep,
            surface,
        }
    }
}

impl Beta for SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
    fn beta(
        &self,
        a_on_d: f64,
        a_on_c: f64,
        _c_on_b: f64,
        _a_on_r: f64,
        _phis: &[f64],
    ) -> Vec<f64> {
        vec![
            self.deep.interp(a_on_d, a_on_c),
            self.surface.interp(a_on_d, a_on_c),
        ]
    }

    fn area(&self, a_on_d: f64, a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        2.0 * a_on_d * a_on_c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_double_corner_crack_hole_newman84() {}
    #[test]
    fn test_centre_crack_tension_fedderson66() {
        let centre = CentreCrackTensionFedderson66 {};

        let phis = vec![0.0f64, FRAC_PI_2];
        assert!(centre.beta(0.0, 0.0, 0.0, 0.0, &phis)[0] - 1.0 < std::f64::EPSILON);
    }

    #[test]
    fn test_centre_crack_tension_koiter65() {
        let centre = CentreCrackTensionKoiter65 {};
        let phis = vec![0.0, FRAC_PI_2];

        assert!(centre.beta(0.0, 0.0, 0.0, 0.0, &phis)[0] - 1.0 < std::f64::EPSILON);
    }

    #[test]
    fn test_mcdonald07() {
        let corner = CornerCrackConstrainedTensionMcdonald07::new();
        let phis = vec![0.0, FRAC_PI_2];

        assert!((corner.beta(0.0, 0.0, 0.0, 0.0, &phis)[0] - 0.709411963).abs() < std::f64::EPSILON);
        assert!((corner.beta(0.00010 / 25.0e-3, 0.0, 0.0, 0.0, &phis)[0] - 0.709411963).abs() < std::f64::EPSILON);
        assert!((corner.beta(0.00987285 / 25.0e-3, 0.0, 0.0, 0.0, &phis)[0] - 1.4485189776434797).abs() < 1e-3);
    }

    #[test]
    // check that the if statement arms a_on_c <> 1 produce the same answer
    fn check_quarter_arms() {
        let offset = 1e-6;

        let quarter = QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {};
        let phis = vec![0.0, FRAC_PI_2];

        let lower = quarter.beta(0.0, 1.0 - offset, 0.0, 0.0, &phis);
        let upper = quarter.beta(0.0, 1.0 + offset, 0.0, 0.0, &phis);

        assert!((lower[0] - upper[0]).abs() < 1e-6);

        let semi = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {};

        let lower = semi.beta(0.0, 1.0 - offset, 0.0, 0.0, &phis);
        let upper = semi.beta(0.0, 1.0 + offset, 0.0, 0.0, &phis);

        assert!((lower[0] - upper[0]).abs() < 1e-6);
    }

    #[test]
    // Python
    // np.sqrt(((2/np.pi) / a_on_b) * np.tan(a_on_b * np.pi/2)) *
    // (0.752 + 2.02 * a_on_b + 0.37 * (1 - np.sin(a_on_b * np.pi/2))**3) / np.cos(a_on_b * np.pi/2))
    fn check_single_edge_notched_tension_tada73() {
        let single = SingleSidedEdgeCrackTensionTada73 {};
        let phis = vec![0.0, FRAC_PI_2];

        println!(
            "tada73 sent(0.0) = {}",
            single.beta(0.0, 0.0, 0.0, 0.0, &phis)[0]
        );

        let single = SingleSidedEdgeCrackTensionTada73 {};

        assert!((single.beta(0.0, 0.0, 0.0, 0.0, &phis)[0] - 1.122).abs() < 1e-6);
    }

    #[test]
    fn check_finite_quarter_beta() {
        let phis = vec![0.0, FRAC_PI_2];
        let quarter = QuarterCircularCornerCrackFinitePlateTensionMurakami87 {};

        assert!((quarter.beta(0.0f64, 0.0, 0.0, 0.0, &phis)[0] - 0.671604).abs() < 1e-6);
        assert!((quarter.beta(0.5f64, 0.0, 0.0, 0.0, &phis)[0] - 0.724126).abs() < 1e-6);
        assert!((quarter.beta(1.0f64, 0.0, 0.0, 0.0, &phis)[0] - 0.7622).abs() < 1e-6);
    }

    #[test]
    fn check_semi_elliptical_infinite_anderson05() {
        let phis = vec![0.0, FRAC_PI_2];
        let semi = SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {};

        assert!((semi.beta(0.0, 1.0, 0.0, 0.0, &phis)[0] - 0.662541348868913).abs() < std::f64::EPSILON);
        assert!((semi.beta(0.0, 0.5, 0.0, 0.0, &phis)[0] - 0.8959634744787476).abs() < std::f64::EPSILON);
    }

    #[test]
    fn check_semi_elliptical_newman84() {
        let phis = vec![0.0, FRAC_PI_2];
        let semi = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {};

        assert!(semi.beta(0.0, 1.0, 0.0, 0.0, &phis)[0] - 0.6625413488689131 < std::f64::EPSILON);
        assert!(semi.beta(0.0, 0.5, 0.0, 0.0, &phis)[0] - 0.8959634744787476 < std::f64::EPSILON);
    }

    #[test]
    fn print_results() {
        let phis = vec![0.0, FRAC_PI_2];
        let semi = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {};
        println!(
            "SemiEllipticalSurfaceCrackFinitePlateTensionNewman84(1.0, [0, pi/2]) {:?}",
            semi.beta(1.0, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SemiEllipticalSurfaceCrackFinitePlateTensionNewman84(0.5, [0, pi/2]) {:?}",
            semi.beta(0.5, 0.0, 0.0, 0.0, &phis)
        );

        // equal a and c betas
        println!(
            "SemiEllipticalSurfaceCrackFinitePlateTensionNewman84(0.82, [0, pi/2]) {:?}",
            semi.beta(0.82, 0.0, 0.0, 0.0, &phis)
        );

        println!(
            "SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05(1.0, [0, pi/2]) {:?}",
            semi.beta(1.0, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05(0.5, [0, pi/2]) {:?}",
            semi.beta(0.5, 0.0, 0.0, 0.0, &phis)
        );

        let quarter = QuarterBroek86 {};
        println!("Quarterbroek86 {:?}", quarter.beta(0.0, 0.0, 0.0, 0.0, &phis));

        let quarter = QuarterCircularCornerCrackFinitePlateTensionMurakami87 {};
        println!(
            "QuarterCircularCornerCrackFinitePlateTensionMurakami87(a_on_t = 0) {:?}",
            quarter.beta(0.0, 0.0, 0.0, 0.0, &phis)
        );
        println!("QuarterEllipticalCornerCrackFinitePlateTensionNewman84(a_on_c=1.0, a_on_d=0.0, c_on_b=0, [0, pi/2]) {:?}", quarter.beta(1.0, 0.0, 0.0, 0.0, &phis));
        println!("QuarterEllipticalCornerCrackFinitePlateTensionNewman84(a_on_c=0.5, a_on_d=0.0, c_on_b=0, [0, pi/2]) {:?}",   quarter.beta(0.5, 0.0, 0.0, 0.0, &phis));

        let elliptical = EllipticalEmbeddedCrackFinitePlateTensionNewman84 {};
        println!("EllipticalEmbeddedCrackFinitePlateTensionNewman84(a_on_c=1.0, a_on_d=0.0, c_on_b=0, [0, pi/2]) {:?}",  elliptical.beta(1.0, 0.0, 0.0, 0.0, &phis));

        let single = SingleSidedEdgeCrackTensionTada73 {};
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 0) {:?}",
            single.beta(0.0, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 1e-6) {:?}",
            single.beta(1e-6, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 0.1) {:?}",
            single.beta(0.1, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 0.5) {:?}",
            single.beta(0.5, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 1.2) {:?}",
            single.beta(1.2, 0.0, 0.0, 0.0, &phis)
        );

        let double = DoubleSidedEdgeCrackTensionTada73 {};
        println!(
            "DoubleSidedEdgeCrackTensionTada73(a_on_d = 0) {:?}",
            double.beta(0.0, 0.0, 0.0, 0.0, &phis)
        );

        let compact_tension_tada73 = CompactCoupon {
            width: 1.0,
            depth: 1.0,
        };
        println!(
            "compact_tension_tada73(a_on_d = 0.5, b=1, w = 1) {:?}",
            compact_tension_tada73.beta(0.5, 1.0, 1.0, 0.0, &phis)
        );
    }

    #[test]
    fn check_semi_elliptical_surface_crack_round_bar_bending_shin04() {
        let phis = vec![0.0, FRAC_PI_2];
        let semi = SemiEllipticalSurfaceCrackRoundBarBendingShin04::new();

        let result = semi.beta(0.067, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.963).abs() < 1e-3);

        let result = semi.beta(0.2, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.872).abs() < 1e-3);

        let result = semi.beta(0.5, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 1.19672).abs() < 1e-3);

        let result = semi.beta(0.8, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 5.140).abs() < 1e-3);

        let result = semi.beta(0.067, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.576).abs() < 1e-3);

        let result = semi.beta(0.2, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.506).abs() < 1e-3);

        let result = semi.beta(0.5, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.524).abs() < 1e-3);

        let result = semi.beta(0.8, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 1.737).abs() < 1e-3);

        let result = semi.beta(0.067, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.81).abs() < 1e-3);

        let result = semi.beta(0.2, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.761).abs() < 1e-3);

        let result = semi.beta(0.5, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 1.0348).abs() < 1e-3);

        let result = semi.beta(0.8, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 4.454).abs() < 1e-3);
    }

    #[test]
    fn manual_check_print_round_bar_table_murakami87() {
        // prints the beta table from murakami P. 658.
        // The first row of the table for a/d = 0 is correct but there is a slight error thereafter.
        // This implies there is an error with the a/d effect (a/d > 0)

        // This is the generated version of the beta table using the code
        // below. The results in this table are slightly different from the
        // that provided on P.658 of Murakami. Since this beta is based on the
        // following three betas:
        //
        // semi_elliptical_surface_crack_round_bar_tension_murakami87(*a_on_d, *crack.a_on_c);
        // edge_crack_strip_bending_murakami87(*a_on_d);
        // edge_crack_strip_tension_murakami87(*a_on_d);
        //
        // then there is either a mistake in my implementation of one or more of
        // these three or in the table that was produced by Murakami.

        // However, the formula given by murakami must be in error since the
        // bending of a crack must start with term of 1.12, but this is for the
        // equation containing the crack aspect ratio beta = b/a which is the
        // aspect ratio of the crack.

        //                                           a/c
        // a/d  | 0.000  0.100  0.200  0.300  0.400  0.500  0.600  0.700  0.800  0.900  1.000
        //      -----------------------------------------------------------------------------
        // 0    | 1.123  1.092  1.048  0.996  0.940  0.884  0.829  0.778  0.733  0.694  0.661
        // 0.01 | 1.114  1.084  1.040  0.989  0.933  0.877  0.823  0.773  0.728  0.689  0.656
        // 0.02 | 1.105  1.074  1.031  0.980  0.925  0.869  0.816  0.766  0.721  0.683  0.650
        // 0.03 | 1.094  1.064  1.021  0.970  0.916  0.861  0.808  0.758  0.714  0.676  0.643
        // 0.04 | 1.082  1.052  1.010  0.960  0.906  0.852  0.799  0.750  0.707  0.669  0.637
        // 0.05 | 1.070  1.041  0.999  0.949  0.896  0.842  0.790  0.742  0.699  0.661  0.630
        // 0.06 | 1.058  1.028  0.987  0.938  0.885  0.832  0.781  0.733  0.691  0.654  0.622
        // 0.07 | 1.045  1.016  0.975  0.927  0.875  0.822  0.771  0.724  0.682  0.646  0.615
        // 0.08 | 1.032  1.003  0.963  0.915  0.864  0.812  0.762  0.715  0.674  0.637  0.607
        // 0.09 | 1.018  0.990  0.950  0.903  0.853  0.801  0.752  0.706  0.665  0.629  0.599
        // 0.1  | 1.005  0.977  0.938  0.892  0.842  0.791  0.742  0.697  0.656  0.621  0.591
        // 0.11 | 0.992  0.965  0.926  0.880  0.831  0.781  0.733  0.688  0.648  0.613  0.584
        // 0.12 | 0.979  0.952  0.914  0.869  0.820  0.771  0.723  0.679  0.640  0.605  0.576
        // 0.13 | 0.967  0.940  0.902  0.858  0.810  0.761  0.714  0.670  0.631  0.597  0.569
        // 0.14 | 0.955  0.928  0.891  0.847  0.799  0.751  0.705  0.662  0.623  0.590  0.562
        // 0.15 | 0.943  0.917  0.880  0.836  0.789  0.742  0.696  0.654  0.616  0.583  0.555
        // 0.16 | 0.931  0.905  0.869  0.826  0.780  0.733  0.688  0.646  0.608  0.575  0.548
        // 0.17 | 0.920  0.895  0.859  0.816  0.770  0.724  0.679  0.638  0.601  0.569  0.541
        // 0.18 | 0.909  0.884  0.849  0.807  0.761  0.716  0.671  0.630  0.594  0.562  0.535
        // 0.19 | 0.899  0.874  0.839  0.797  0.753  0.707  0.664  0.623  0.587  0.556  0.529
        // 0.2  | 0.889  0.864  0.830  0.789  0.744  0.700  0.656  0.616  0.580  0.549  0.523
        // 0.21 | 0.879  0.855  0.821  0.780  0.736  0.692  0.649  0.610  0.574  0.543  0.517
        // 0.22 | 0.870  0.846  0.812  0.772  0.728  0.685  0.642  0.603  0.568  0.538  0.512
        // 0.23 | 0.861  0.837  0.804  0.764  0.721  0.677  0.636  0.597  0.562  0.532  0.506
        // 0.24 | 0.852  0.829  0.795  0.756  0.714  0.671  0.629  0.591  0.556  0.527  0.501
        // 0.25 | 0.844  0.820  0.787  0.748  0.706  0.664  0.623  0.585  0.551  0.521  0.496
        
        let table_beta = SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {};
        println!("{}", table_beta.as_table());
    }

    #[test]
    fn manual_check_print_round_bar_table_murakami86() {
        // prints the beta table from Murakami and Tsuru 86
        let table_beta = SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {};
        println!("{}", table_beta.as_table());
    }

    #[test]
    // check that for the same conditions the quarter elliptical
    // always produces a higher beta than the semi elliptical
    fn check_semi_quarter_elliptical_crack_newman84() {
        let phis = vec![0.0, FRAC_PI_2];

        let a_on_ds = [0.0f64, 0.2, 0.5, 0.8];
        let a_on_cs = [0.2f64, 0.5, 0.8];
        let phis = [0.0];
        let c_on_b = 0.0;
        let a_on_r = 0.0;

        for &a_on_d in &a_on_ds {
            for &a_on_c in &a_on_cs {
                println!("a_on_d {}, a_on_c {}", a_on_d, a_on_c);

                let semi = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {};
                let s = semi.beta(a_on_d, a_on_c, c_on_b, a_on_r, &phis);

                let quarter = QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {};
                let q = quarter.beta(a_on_d, a_on_c, c_on_b, a_on_r, &phis);

                println!("semi {:?}, quarter {:?}", s, q);
                assert!(q > s);
            }
        }
    }
}
