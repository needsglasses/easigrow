//! Calculations for a compact tension coupon.
//!
//! Equations are take from the reference:
//! Back-face strain compliance relation for compact specimens for wide
//! range in crack lengths
//! J.C. Newman Jr., Y. Yamada, M.A. James
//! Engineering Fracture Mechanics 78 (2011) 2707–2711

#[allow(dead_code)]
/// Calculate the crack length from backface strain measurements.
fn cracklength_from_backface_strain(strain: f64, youngs_modulus: f64, depth: f64, load: f64) -> f64{
    let a = youngs_modulus * strain * depth / load;
    let u = 1.0 /(a.sqrt() + 1.0);

    let a0 = 1.0033;
    let a1 = -2.35;
    let a2 = 1.3694;
    let a3 = -15.294;
    let a4 = 63.182;
    let a5 = -74.42;

    let a_on_d =  a0 + a1 * u + a2 * u.powi(2) + a3 * u.powi(3) + a4 * u.powi(4) + a5 * u.powi(5);

    a_on_d * depth
}


/// Calculate the back face strain knowing the crack size.
pub fn backface_strain(a_on_d: f64, load: f64, youngs_modulus: f64, depth:f64, width:f64) -> f64 { 
    let compliance = (1.41 - 1.462 * a_on_d + 20.45 * a_on_d.powi(2) 
                      - 26.83 * a_on_d.powi(3) + 11.45 * a_on_d.powi(4)) / (1.0 - a_on_d).powi(2);

    let strain = compliance * load /(youngs_modulus * depth * width);

    strain
}


/// Displacement of a front face clip gauge.
///
/// from ASTM 647
pub fn clipgauge_displacement(a_on_d: f64, load: f64, youngs_modulus: f64, depth:f64, width:f64) -> f64 { 
    // this equation needs to be inverted to solve for u
    let v = 0.0; // clipgauge displacement
    let u = ((youngs_modulus * v * width / load).sqrt() + 1).recip();
       
    let a_on_d = 1.0010 - 4.6695 * u + 18.460 * u.powi(2) 
        - 236.82 * u.powi(3) + 1214.9 * u.powi(4) - 2143.6 * u.powi(5);
    
    a_on_d
}

#[test]
/// From figure 2 in the paper a ratio of a/w = 0.45 has a compliance of 10.
/// compliance is E(epsilon W) B/P
fn check_compliance_backface() {
    let strain = backface_strain(1.0, 1.0, 1.0, 1.0, 1.0);
    
    println!("compliance check {}", strain - 10.0);
}
