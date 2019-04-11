//! Collection of da/dn equations for calculating crack growth.
//!
//! So far these subroutines are history independant and only depend
//! on $\Delta K$ and $R$.

use table;
use std::f64::consts::FRAC_PI_2;
use std::{fmt, process, f64};
use log::{info, warn, error};

/// Defines the state of the crack for any dadn equation that requires
/// some sort of memory.
pub struct CrackState {
    // length of the crack
    pub a: f64
}

pub enum DadnEqn {
    Nasgro,
    Forman,
    Paris,
    Walker,
    Burchill,
    Hartman,
    White,
    Kujawski,
}

/// data for White equation
/// $$ dadn =
#[derive(Debug, Clone)]
pub struct White {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    kic: f64,
    cite: String,
    units: String,
}

/// data for Kujawski equation
/// $$ dadn =
#[derive(Debug, Clone)]
pub struct Kujawski {
    c: f64,
    m: f64,
    alpha: f64,
    cite: String,
    units: String,
}

// Burchill equation
// dadn = C Kmax^m - D Kmin^n
#[derive(Debug, Clone)]
pub struct Burchill {
    c: f64,
    m: f64,
    d: f64,
    n: f64,
    cite: String,
    units: String,
}

/// data for Paris equation
/// $$ dadn = C \Delta K^m $$
#[derive(Debug, Clone)]
pub struct Paris {
    c: f64,
    m: f64,
    cite: String,
    units: String,
}

/// Forman equation
/// $$ da/dn = C \Delta K^n / ((1-R) * `K_f`  - \Delta K) $$
#[derive(Debug, Clone)]
pub struct Forman {
    c: f64,
    n: f64,
    kf: f64,
    cite: String,
    units: String,
}

/// Walker Equation
/// $$ dadn = C \Delta K^m  $$
#[derive(Debug, Clone)]
pub struct Walker {
    c: f64,
    m: f64,
    n: f64,
    cite: String,
    units: String,
}

/// Hartman-Schijve Variant (Jones-Molent)
/// The influence of cyclic stress intensity threshold on fatigue life scatter
/// L. Molent and R. Jones
/// International Journal of Fatigue
///Volume 82, Part 3, January 2016, Pages 748–756
#[derive(Debug, Clone)]
pub struct Hartman {
    d: f64,
    k_thr: f64,
    a: f64,
    alpha: f64,
    cite: String,
    units: String,
}

#[derive(Debug, Clone)]
pub struct Nasgro {
    // ratio of maximum far field stress to flow stress sigma0
    smax_on_sigma0: f64,
    // constraint factor
    alpha: f64,
    // fracture toughness
    k_crit: f64,
    // delta k threshold at R=0
    deltak0: f64,
    // curve control coefficient for different values of R, equals 0 for negative R, equals 1 for R>=0.
    cth: f64,
    cth_minus: f64,
    p: f64,
    // emperical crack coefficient
    q: f64,
    // emperical crack coefficient
    c: f64,
    // emperical crack coefficient
    n: f64,
    // intrinsic crack size (m)
    a_intr: f64,
    cite: String,
    units: String,
}

/// Trait for `DaDn` function that can be used for crack growth
pub trait DaDn {
    /// Display the crack growth equation
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result;

    /// Update the equation with new variables
    fn update(&mut self, x: &[f64]);

    /// Extract the variables from the equation and return as a simple vector.
    fn variables(&self) -> Vec<f64>;

    /// dadn equation gives the instaneous rate of crack growth
    ///
    /// Currently we also pass it the crack length which is currently
    /// not used by all equations but it could be generalise to be a
    /// general memory state parameter.
    fn dadn(&self, kmin: f64, kmax: f64, state: CrackState) -> f64;
}

pub struct Closure<T: ?Sized> {
    a_start: f64,
    a_end: f64,
    closure_start: f64,
    closure_end: f64,
    eqn: Box<T>,
}

impl<T: ?Sized> Closure<T> {
    // Simple linear model of closure.
    pub fn new(
        model: Box<T>,
        a_start: f64,
        a_end: f64,
        closure_start: f64,
        closure_end: f64,
    ) -> Closure<T> {
        Closure {
            a_start,
            a_end,
            closure_start,
            closure_end,
            eqn: model,
        }
    }

    fn closure(&self, a: f64) -> f64 {
        self.closure_start
            + (a - self.a_start) * (self.closure_start - self.closure_end)
                / (self.a_end - self.a_start)
    }
}

/// calculates a new kmin based on a simple closure model
impl<T: DaDn + ?Sized> DaDn for Closure<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "")
    }

    fn update(&mut self, x: &[f64]) {
        self.eqn.update(x);
        let n = x.len();

        self.closure_start = x[n - 2];
        self.closure_end = x[n - 1];
    }

    fn variables(&self) -> Vec<f64> {
        vec![]
    }

    fn dadn(&self, kmin: f64, kmax: f64, state: CrackState) -> f64 {
        let kmin_with_closure = (kmax * self.closure(state.a)).max(kmin);

        self.eqn.dadn(kmin_with_closure, kmax, state)
    }
}

impl fmt::Display for DaDn {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)
    }
}

/// White equation.
/// An extension of the Forman equation to create a dadn curve made from
/// a combination of a cubic with an asymptotic curve for kic effect.
/// Used in the specification for MSMP3.

impl White {
    pub fn new(x: &[f64], cite: String) -> White {
        White {
            a: x[0],
            b: x[1],
            c: x[2],
            d: x[3],
            e: x[4],
            f: x[5],
            kic: x[6],
            cite,
            units: String::from("m"),
        }
    }
}

impl DaDn for White {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:#?}", self);
        write!(f, 
"#  da/dN ({units}) = exp[({a:e} * ΔKeff^3 - {b:e} * ΔKeff^2 + {c:e} * ΔKeff - {d:e}) 
#              + (dkic - ΔK)^{f}] [White14:{cite}]
#  where dkic = {kic:e} (1 - R) and ΔKeff = ΔK / (1 - R)^{e}",
               units=self.units, a=self.a, b=self.b, c=self.c, d=self.d, e=self.e, kic=self.kic, f=-self.f, cite=self.cite )
    }

    fn update(&mut self, x: &[f64]) {
        self.a = x[0];
        self.b = x[1];
        self.c = x[2];
        self.d = x[3];
        self.e = x[4];
        self.f = x[5];
        self.kic = x[6];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.a, self.b, self.c, self.d, self.e, self.f, self.kic]
    }

    // note the signs of the b,d and f variables has been change to make all the coefficients positive
    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        // this function does not work well if the rmax variable is
        // included in the optimisation so we set it directly
        // let rmax = coeffs[7];
        let rmax = 0.8; // hardwired

        let delta_k = kmax - kmin.max(0.0);
        let r = (kmin.max(0.0) / kmax).min(rmax);

        let dkic = self.kic * (1.0 - r);
        let d_ke = (delta_k / (1.0 - r).powf(self.e)).ln();

        // simple cubic curve
        let main = self.a * d_ke.powi(3) - self.b * d_ke.powi(2) + self.c * d_ke - self.d;

        // extra bit due to approaching fracture toughness
        let bit = (dkic - delta_k).max(0.01).powf(-self.f);

        if delta_k > 0.0 {
            (main + bit).exp()
        } else {
            0.0
        }
    }
}

/// Kujawski equation
/// parameters = [c, n, k]
impl Kujawski {
    pub fn new(x: &[f64], cite: String) -> Kujawski {
        Kujawski {
            c: x[0],
            m: x[1],
            alpha: x[2],
            cite,
            units: String::from("m"),
        }
    }
}

impl DaDn for Kujawski {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * [Kmax ^ {alpha} * ΔK^ {{1 - {alpha}}}]^{m}  [Kujawski01]",
            units = self.units,
            c = self.c,
            m = self.m,
            alpha = self.alpha,
        )
    }

    fn update(&mut self, x: &[f64]) {
        self.c = x[0];
        self.m = x[1];
        self.alpha = x[2];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.c, self.m, self.alpha]
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        let delta_k = kmax.max(0.0) - kmin.max(0.0);
        self.c * (kmax.max(0.0).powf(self.alpha) * delta_k.powf(1.0 - self.alpha)).powf(self.m)
    }
}

/// Burchill equation
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
impl Burchill {
    pub fn new(x: &[f64], cite: String) -> Burchill {
        Burchill {
            c: x[0],
            m: x[1],
            d: x[2],
            n: x[3],
            cite,
            units: String::from("m"),
        }
    }
}

impl DaDn for Burchill {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * Kmax ^ {m} - {d:e} Kmin ^ {n}]  [Burchill17]",
            units = self.units,
            c = self.c,
            m = self.m,
            d = self.d,
            n = self.n
        )
    }

    fn update(&mut self, x: &[f64]) {
        self.c = x[0];
        self.m = x[1];
        self.d = x[2];
        self.n = x[3];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.c, self.m, self.d, self.n]
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        (self.c * kmax.max(0.0).powf(self.m) - self.d * kmin.max(0.0).powf(self.n))
    }
}

/// Nasgro  equation
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
/// Ref: AFGROW users guide and technical manual
/// James A. Harter
/// AFRL-VA-WP-TR-1999-3016
/// Feb 1999
impl Nasgro {
    pub fn new(x: &[f64], cite: String) -> Nasgro {
        Nasgro {
            smax_on_sigma0: x[0],
            alpha: x[1],
            k_crit: x[2],
            deltak0: x[3],
            cth: x[4],
            cth_minus: x[5],
            p: x[6],
            q: x[7],
            c: x[8],
            n: x[9],
            a_intr: x[10],
            cite,
            units: String::from("m"),
        }
    }
}

// According to the AFGROW manual, the variable
// Smax_on_simga is the ratio of the maximum applied (far field) stress
// to the flow stress (of the material).  This does not make a lot of
// sense as they really have nothing to do with each other.
impl DaDn for Nasgro {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        _e = writeln!(
            f,
            "#  da/dN ({units}) = {c:e} * (((1 - f) / (1 - R)) ΔK ^ {n}) * G  [Nasgro99]",
            units = self.units,
            c = self.c,
            n = self.n
        );
        writeln!(
            f,
            "#  G = (1 - (ΔK_th / ΔK)) ^ {p} /  (1 - (Kmax / Kcrit)) ^ {q}   ",
            p = self.p,
            q = self.q
        )
    }

    fn update(&mut self, x: &[f64]) {
        self.smax_on_sigma0 = x[0];
        self.alpha = x[1];
        self.k_crit = x[2];
        self.deltak0 = x[3];
        self.cth = x[4];
        self.cth_minus = x[5];
        self.p = x[6];
        self.q = x[7];
        self.c = x[8];
        self.n = x[9];
        self.a_intr = x[10];
    }

    fn variables(&self) -> Vec<f64> {
        vec![
            self.smax_on_sigma0,
            self.alpha,
            self.k_crit,
            self.deltak0,
            self.cth,
            self.cth_minus,
            self.p,
            self.q,
            self.c,
            self.n,
            self.a_intr,
        ]
    }

    fn dadn(&self, kmin: f64, kmax: f64, state: CrackState) -> f64 {
        let r = kmin / kmax;
        let delta_k = kmax - kmin;

        // closure constants
        let a0 = (0.825 - 0.34 * self.alpha + 0.05 * self.alpha.powi(2))
            * (FRAC_PI_2 * self.smax_on_sigma0)
                .cos()
                .powf(1.0 / self.alpha);
        let a1 = (0.415 - 0.071 * self.alpha) * self.smax_on_sigma0;
        let a3 = 2.0 * a0 + a1 - 1.0;
        let a2 = 1.0 - a0 - a1 - a3;

        // closure level f
        let f = if r >= 0.0 {
            (a0 + a1 * r + a2 * r.powi(2) + a3 * r.powi(3)).max(r)
        } else if (-2.0 <= r) && (r < 0.0) {
            a0 + a1 * r
        } else {
            a0 - 2.0 * a1
        };

        info!("Nasgro: {} {} {} {} {}", a0, a1, a3, a2, f);
        // a_int = intrinsic crack size, typically 38.1e-6 (m)
        let deltak_th = self.deltak0 * (state.a / (state.a + self.a_intr)).sqrt()
            / ((1.0 - f) / ((1.0 - a0) * (1.0 - r))).powf(1.0 + self.cth * r);

        info!("nasgro: deltak_th {}", deltak_th);
        let num = (1.0 - (deltak_th.min(delta_k - 1e-6) / delta_k)).powf(self.p);
        let denom = (1.0 - (kmax / self.k_crit)).powf(self.q);
        let dadn = self.c * (((1.0 - f) / (1.0 - r)) * delta_k).powf(self.n) * num / denom;
        info!("nasgro: dadn {} {} {}", dadn, num, denom);
        dadn
    }
}

/// Forman equation
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
impl Forman {
    pub fn new(x: &[f64], cite: String) -> Forman {
        Forman {
            c: x[0],
            n: x[1],
            kf: x[2],
            cite,
            units: String::from("m"),
        }
    }
}

impl DaDn for Forman {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * ΔK ^ {n} / [(1 - R) {kf} - ΔK]  [Forman67]",
            units = self.units,
            c = self.c,
            n = self.n,
            kf = self.kf
        )
    }

    fn update(&mut self, x: &[f64]) {
        self.c = x[0];
        self.n = x[1];
        self.kf = x[2];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.c, self.n, self.kf]
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        let delta_k = kmax - kmin.max(0.0);
        let r = kmin.max(0.0) / kmax;

        if delta_k > 0.0 {
            self.c * delta_k.powf(self.n) / ((1.0 - r) * self.kf - delta_k)
        } else {
            0.0
        }
    }
}

/// Paris equation
/// da/dN = c (\Delta K)^m
/// parameters = [c, m]
impl Paris {
    pub fn new(x: &[f64], cite: String) -> Paris {
        Paris {
            c: x[0],
            m: x[1],
            cite,
            units: String::from("m"),
        }
    }
}

impl DaDn for Paris {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({}) = {:e} * ΔK ^ {} [Paris63]",
            self.units, self.c, self.m
        )
    }

    fn update(&mut self, x: &[f64]) {
        self.c = x[0];
        self.m = x[1];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.c, self.m]
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        let delta_k = kmax.max(0.0) - kmin.max(0.0);

        self.c * delta_k.powf(self.m)
    }
}

/// Hartman-Schijve equation
/// c [(\Delta K - kth) / \sqrt{1 - (kmax/a)}]^n
/// parameters = [c, n, kth, a]
impl Hartman {
    pub fn new(x: &[f64], cite: String) -> Hartman {
        Hartman {
            d: x[0],
            k_thr: x[1],
            a: x[2],
            alpha: x[3],
            cite,
            units: String::from("m"),
        }
    }
}

impl DaDn for Hartman {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(f, "#  da/dN ({units}) =  {d:e} [(ΔK - {k_thr}) / sqrt (1 - kmax/{a})]^{alpha} [Hartman70]",
               units=self.units, d=self.d, k_thr=self.k_thr, a=self.a, alpha=self.alpha)
    }

    fn update(&mut self, x: &[f64]) {
        self.d = x[0];
        self.k_thr = x[1];
        self.a = x[2];
        self.alpha = x[3];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.d, self.k_thr, self.a, self.alpha]
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        let delta_k = kmax - kmin.max(0.0);
        let r = kmin.max(0.0) / kmax;

        let kmax = delta_k / (1.0 - r);

        if kmax > self.a {
            warn!("***Warning: A kmax of {} is > {} and therefore cannot be square-rooted in the Hartman-Schijve equation", kmax, self.a);
        }

        if delta_k > self.k_thr {
            self.d
                * ((delta_k - self.k_thr).max(0.0) / (1.0 - (kmax / self.a)).max(0.0).sqrt())
                    .powf(self.alpha)
        } else {
            0.0
        }
    }
}

/// Walker equation
///
/// da/dn = c * (\Delta K * (1 - R)^(m - 1))^n
/// parameters = [C, n, m]
impl Walker {
    pub fn new(x: &[f64], cite: String) -> Walker {
        Walker {
            c: x[0],
            m: x[1],
            n: x[2],
            cite,
            units: String::from("m"),
        }
    }
}

/// In this implementation R is not limited to a minimium of 0.0 .
impl DaDn for Walker {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * [ΔK / (1-R)^(1 - {m})] ^ {n} [Walker70]",
            units = self.units,
            c = self.c,
            m = self.m,
            n = self.n
        )
    }

    fn update(&mut self, x: &[f64]) {
        self.c = x[0];
        self.m = x[1];
        self.n = x[2];
    }

    fn variables(&self) -> Vec<f64> {
        vec![self.c, self.m, self.n]
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        let delta_k = (kmax - kmin.max(0.0)).max(0.0);
        let r = (kmin / kmax).max(0.0);

        self.c * (delta_k / (1.0 - r).powf(1.0 - self.m)).powf(self.n)
    }
}

impl DaDn for table::Table {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = write!(
            f,
            "  Spline interpolated tabular lookup
#       da/dn:"
        );

        for col in &self.columns {
            let _ = write!(f, " {:12}", col);
        }
        let _ = writeln!(f);
        for i in 0..self.row.len() {
            let _ = write!(f, "#  {:10.3e}: ", (10.0f64).powf(self.row[i]));
            for j in 0..self.values.len() {
                let _ = write!(f, "  {:10} ", self.values[j][i]);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }

    /// copy the optimisation parameters into the delta k values of the table
    fn update(&mut self, params: &[f64]) {
        let mut i = 0;
        for col in self.values.iter_mut() {
            for cell in col.iter_mut() {
                *cell = params[i];
                i += 1;
            }
        }
    }

    // extract the dk columns out and return as a vector
    fn variables(&self) -> Vec<f64> {
        self.values
            .iter()
            .flat_map(|v| v.clone())
            .collect::<Vec<f64>>()
    }
    
    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        info!("kmin {} kmax {}", kmin, kmax);
        let rmax = self.columns[self.columns.len() - 1];
        let rmin = self.columns[0];
        let delta_k = kmax - kmin.max(0.0);
        let r = (kmin / kmax).max(rmin).min(rmax);

        info!("delta k {dk} r limited {r}", dk=delta_k, r=r);
        let interp = self.interp(delta_k, r);
        info!("Interp {}", interp);
        let value = (10.0f64).powf(interp);
        info!("value {}", value);
        value
    }
}

impl DaDn for table::PairTable {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = write!(
            f,
            "  Spline interpolated tabular lookup
#       da/dn:"
        );

        for col in &self.columns {
            let _ = write!(f, " {:12}", col);
        }
        let _ = writeln!(f);
        for i in 0..self.rows.len() {
            for j in 0..self.values.len() {
                let _ = write!(f, "  {:10.3e}: ", (10.0f64).powf(self.rows[j][i]));
                let _ = write!(f, "  {:10} ", self.values[j][i]);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }

    /// copy the optimisation parameters into the delta k values of the table
    fn update(&mut self, params: &[f64]) {
        let mut i = 0;
        for col in self.values.iter_mut() {
            for cell in col.iter_mut() {
                *cell = params[i];
                i += 1;
            }
        }
    }

    // extract the dk columns out and return as a vector
    fn variables(&self) -> Vec<f64> {
        self.values
            .iter()
            .flat_map(|v| v.clone())
            .collect::<Vec<f64>>()
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        info!("kmin {} kmax {}", kmin, kmax);
        let rmax = self.columns[self.columns.len() - 1];
        let rmin = self.columns[0];
        let delta_k = kmax - kmin.max(0.0);
        let r = (kmin / kmax).max(rmin).min(rmax);

        info!("delta k {dk} r limited {r}", dk=delta_k, r=r);
        let interp = self.interp(delta_k, r);
        info!("Interp {}", interp);
        let value = (10.0f64).powf(interp);
        info!("value {}", value);
        value
    }
}

/// Create a function closure for a dadn model.
///
/// The dadn model name is specified by two parts separated by a colon
/// e.g.  model:material e.g. forman:2024t3-sheet.
/// The data component must be one of the pre-existing models. But if
/// parameters are supplied these will be installed into the dadn
/// model and used.
pub fn make_model(model_name: &str, params: &[f64], cite: String) -> Box<DaDn> {
    let components: Vec<&str> = model_name.split(':').collect();
    let dadn = components[0];

    if dadn.contains("white") {
        Box::new(White::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("forman") {
        Box::new(Forman::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("paris") {
        Box::new(Paris::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("hartman") {
        Box::new(Hartman::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("walker") {
        Box::new(Walker::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("burchill") {
        Box::new(Burchill::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("kujawski") {
        Box::new(Kujawski::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("nasgro") {
        Box::new(Nasgro::new(&params, cite)) as Box<DaDn>
    } else if dadn.contains("file") {
        let table = table::Table::read_file(components[1], true);
        // put the params back into the table
        let mut log_table = table::Table::new(
            table.columns,
            table.row.iter().map(|x| x.log10()).collect::<Vec<f64>>(),
            table.values,
            true,
        );
        if !params.is_empty() {
            if params.len() == log_table.values.iter().fold(0, |sum, x| sum + x.len()) {
                log_table.update(&params);
            } else {
                error!(
                    "Error: the number of parameters {} != the no. of table values {}",
                    params.len(),
                    log_table.values.len()
                );
                process::exit(1);
            }
        }

        Box::new(log_table) as Box<DaDn>
    } else {
        error!("Error: Unknown dadn equation: {:?}", model_name);
        process::exit(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use table;
    use material;

    #[test]
    fn check_kujawski() {
        let materials = material::get_all_dadns();
        let kujawski_eqn = match materials.iter().find(|m| m.name == "kujawski:default") {
            Some(m) => &m.eqn,
            None => panic!(),
        };

        // 6 ksi sqrt(in) gives dadn = 1.25e-6 in
        let da = kujawski_eqn.dadn(0.0, 6.6, CrackState {a: 0.0});
        println!("kujawski da/dn {}", da);
        assert!((da - 2.87496e-8).abs() < 1.0e-10);

        let kujawski_eqn = Kujawski::new(&[1e-10, 3.0, 0.25], String::from("test"));
        let da = kujawski_eqn.dadn(0.0, 6.6, CrackState {a: 0.0});
        println!("kujawski da/dn {}", da);
        assert!((da - 2.87496e-8).abs() < 1.0e-10);

        let kujawski_eqn = Kujawski::new(&[1e-10, 3.0, 0.25], String::from("test"));
        let da = kujawski_eqn.dadn(3.0, 6.6, CrackState {a: 0.0});
        println!("kujawski da/dn {}", da);
        assert!((da - 0.735_086_7e-8).abs() < 1.0e-10);
       
    }

    #[test]
    fn check_tabular() {
        let table = table::Table::new(
            vec![0.0, 1.0],
            // already logged.
            vec![-8.0, -7.0, -6.0, -5.0],
            // these are the columns
            vec![vec![10.0, 15.0, 20.0, 25.0], vec![5.0, 7.5, 10.0, 12.5]],
            true,
        );

        assert!((table.dadn(0.0, 10.0, CrackState {a: 0.0}) - 1e-8).abs() < std::f64::EPSILON);
        println!("dadn: min -5.0, max 10.0");
        assert!((table.dadn(-5.0, 10.0, CrackState {a: 0.0}) - 1e-8).abs() < std::f64::EPSILON);

        println!("dadn: min 0.0, max 12.5");
        assert!((table.dadn(0.0, 12.5, CrackState {a: 0.0}) - 3.162_277_660_168_379e-08).abs() < std::f64::EPSILON);
        assert!((table.dadn(12.5, 25.0, CrackState {a: 0.0}) - 0.000_000_562_341_325_190_349).abs() < std::f64::EPSILON);
    }

    #[test]
    fn check_walker_rratio() {
        // check that changing the r exponent makes no difference when the r ratio is 0.0
        // this makes sure we know which exponent is which
        let w1 = Walker::new(&[1e-8, 0.5, 2.0], String::from("testing"));
        let w2 = Walker::new(&[1e-8, 2.5, 2.0], String::from("testing"));

        let a1 = w1.dadn(0.0, 10.0, CrackState {a: 0.1});
        let a2 = w2.dadn(0.0, 10.0, CrackState {a: 0.1});

        assert!((a1 - a2).abs() < 1e-8);
    }

    #[test]
    fn check_nasgro() {
        let materials = material::get_all_dadns();
        let nasgro_eqn = match materials.iter().find(|m| m.name == "nasgro:default") {
            Some(m) => &m.eqn,
            None => panic!(),
        };

        // 6 ksi sqrt(in) gives dadn = 1.25e-6 in
        let da = nasgro_eqn.dadn(0.0, 6.6, CrackState {a: 0.0});
        println!("nasgro da/dn {}", da);
        assert!(((da - 3.2668e-8).abs() / da) < 1.0e-3);
    }
}
