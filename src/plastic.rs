//! Calculate various plastic zone properties.

use std::f64::consts::PI;
use std::ops::{Add, Sub};

#[derive(Debug, Clone)]
pub struct ZoneWidth {
    pub plane_stress: f64,
    pub plane_strain: f64,
}

impl Add<f64> for ZoneWidth {
    type Output = ZoneWidth;

    fn add(self, offset: f64) -> ZoneWidth {
        ZoneWidth {
            plane_stress: self.plane_stress + offset,
            plane_strain: self.plane_strain + offset,
        }
    }
}

impl Sub<f64> for ZoneWidth {
    type Output = ZoneWidth;

    fn sub(self, offset: f64) -> ZoneWidth {
        ZoneWidth {
            plane_stress: self.plane_stress - offset,
            plane_strain: self.plane_strain - offset,
        }
    }
}

/// Calculate the size of the plastic zone for (plane stress, plane strain)
/// Ref. Anderson P. 485
/// Ref. Tada 1973 P. 1.17
pub fn zone_size(kmax: f64, sigma_yield: f64) -> ZoneWidth {
    // alpha is a constraint factor
    let alpha_plane_stress = 2.0;
    let alpha_plane_strain = 6.0;

    // radius of onset of yield with adjustment of the shape of the yield zone
    let radius_yield = |alpha: f64| (1.0 / (alpha * PI)) * (kmax / sigma_yield).powi(2);

    // the plastic zone width is approximately 2 * radius_yield
    ZoneWidth {
        plane_stress: 2.0 * radius_yield(alpha_plane_stress),
        plane_strain: 2.0 * radius_yield(alpha_plane_strain),
    }
}
