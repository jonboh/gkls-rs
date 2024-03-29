#![warn(
    missing_debug_implementations,
    missing_docs,
    trivial_casts,
    trivial_numeric_casts,
    unused_extern_crates,
    unused_import_braces,
    unused_qualifications,
    unused_results
)]
//! `gkls-rs` is a pure Rust implementation of the algorithm described in
//! [Software for Generation of Classes of Test Functions with Known Local and Global Minima for Global Optimization](https://arxiv.org/abs/1103.2695).

use std::f64::consts::PI;

use itertools::izip;
use thiserror::Error;

#[cfg(feature = "test_cbinding")]
pub(crate) mod c_gkls;

pub(crate) mod ranf;

#[derive(Clone, Debug)]
pub struct Options {
    /// Penalty value of the generated function if x is not in D
    pub max_value: f64,
    /// Value of the machine zero in the floating-point arithmetic
    pub precision: f64,
    /// Default value of the paraboloid minimum
    pub paraboloid_min: f64,
    pub delta_max_value: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_value: 1e+100,
            precision: 1e-10,
            paraboloid_min: 0.0,
            delta_max_value: 10.0,
        }
    }
}

#[derive(Debug, Error)]
pub enum GKLSError {
    #[error("Input had incorrect dimension")]
    Dim,
    #[error("Invalid number of minima")]
    NumMinima,
    #[error("Function number outside range")]
    FuncNumber,
    #[error("Evaluation point outside problem boundary?")]
    Boundary,
    #[error("Invalid GlobalMinValue")]
    GlobalMinValue,
    #[error("Invalid GlobalDist")]
    GlobalDist,
    #[error("Invalid GlobalRadius")]
    GlobalRadius,
    #[error("Out of memory")]
    Memory,
    #[error("Error in derivative evaluation")]
    DerivEval,
    #[error("FloatingPoint error")]
    FloatingPoint,
}

#[derive(Debug, Clone)]
struct Minima {
    local_min: Vec<Vec<f64>>,
    f: Vec<f64>,
    w_rho: Vec<f64>,
    peak: Vec<f64>,
    rho: Vec<f64>,
}

#[derive(Debug, Clone)]
struct GlobalMinima {
    num_global_minima: usize,
    gm_index: Vec<usize>,
}

enum CoincidenceCondition {
    ParabolaMinCoincidence,
    LocalMinCoincidence,
    Ok,
}

/// Defines the domain of the optimization problem
#[derive(Debug, Clone)]
pub struct Domain {
    /// Left boundary vector of D
    left: Vec<f64>,
    /// Right boundary vector of D
    right: Vec<f64>,
}

impl Domain {
    fn default(dim: usize) -> Self {
        Self {
            left: (0..dim).map(|_| -1.0).collect(),
            right: (0..dim).map(|_| 1.0).collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Problem {
    options: Options,
    domain: Domain,
    dim: usize,
    num_minima: usize,
    _global_dist: f64,
    _global_radius: f64,
    _global_value: f64,
    delta: f64,
    minima: Minima,
    _glob: GlobalMinima,
}

impl Default for Problem {
    fn default() -> Self {
        let nf = 1;
        let dim = 2;
        let num_minima = 10;
        let global_value = -1.0;
        let domain = Domain::default(dim);
        let global_dist = std::iter::zip(domain.left.iter(), domain.right.iter())
            .map(|(left, right)| right - left)
            .fold(None, |min, val| match min {
                None => Some(val),
                Some(m) => Some(f64::min(m, val)),
            })
            .unwrap()
            / 3.0; // TODO: beautify
        let global_radius = 0.5 * global_dist;
        Self::new(
            nf,
            Options::default(),
            dim,
            num_minima,
            global_value,
            global_radius,
            global_dist,
        )
        .unwrap()
    }
}

fn norm(x1: &[f64], x2: &[f64]) -> f64 {
    let mut norm_ = 0.0;
    for i in 0..x1.len() {
        norm_ += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    norm_.sqrt()
}

fn coincidence_check(minima: &Minima, num_minima: usize, precision: f64) -> CoincidenceCondition {
    for i in 2..num_minima {
        if norm(&minima.local_min[i], &minima.local_min[0]) < precision {
            return CoincidenceCondition::ParabolaMinCoincidence;
        }
    }
    for i in 1..num_minima - 1 {
        for j in i + 1..num_minima {
            if norm(&minima.local_min[i], &minima.local_min[j]) < precision {
                return CoincidenceCondition::LocalMinCoincidence;
            }
        }
    }
    CoincidenceCondition::Ok
}

impl Problem {
    fn parameters_check(
        dim: usize,
        num_minima: usize,
        domain: &Domain,
        options: &Options,
        global_value: f64,
        global_radius: f64,
        global_dist: f64,
    ) -> Result<(), GKLSError> {
        if dim <= 1 {
            return Err(GKLSError::Dim);
        }
        if num_minima <= 1 {
            return Err(GKLSError::NumMinima);
        }
        for i in 0..dim {
            if domain.left[i] >= domain.right[i] - options.precision {
                return Err(GKLSError::Boundary);
            }
        }
        if global_value >= options.paraboloid_min - options.precision {
            return Err(GKLSError::GlobalMinValue);
        }
        let min_side = std::iter::zip(domain.left.iter(), domain.right.iter())
            .map(|(left, right)| right - left)
            .fold(None, |min, val| match min {
                None => Some(val),
                Some(m) => Some(f64::min(m, val)),
            })
            .unwrap();
        if global_dist >= 0.5f64.mul_add(min_side, -options.precision)
            || global_dist <= options.precision
        {
            return Err(GKLSError::GlobalDist);
        }
        if global_radius >= 0.5f64.mul_add(global_dist, options.precision)
            || global_radius <= options.precision
        {
            return Err(GKLSError::GlobalRadius);
        }
        Ok(())
    }

    fn set_basins(
        num_minima: usize,
        minima: &mut Minima,
        glob: &mut GlobalMinima,
        global_value: f64,
        global_radius: f64,
        options: &Options,
        rng: ranf::Ranf,
    ) -> Result<(), GKLSError> {
        let mut rng = rng;
        let mut temp_min: f64;
        let mut temp_d1: f64;
        let mut temp_d2: f64;
        let mut dist: f64;
        // Set the radii rho(i) of the attraction regions
        for i in 0..num_minima {
            temp_min = f64::MAX;
            for j in 0..num_minima {
                if i != j {
                    let temp_d1 = norm(&minima.local_min[i], &minima.local_min[j]);
                    if temp_d1 < temp_min {
                        temp_min = temp_d1;
                    }
                }
            }
            dist = temp_min / 2.0;
            minima.rho[i] = dist;
        }
        // Adjust the radii of the attraction regions
        minima.rho[1] = global_radius;
        for i in 2..num_minima {
            dist = norm(&minima.local_min[i], &minima.local_min[1])
                - global_radius
                - options.precision;
            if dist < minima.rho[i] {
                minima.rho[i] = dist;
            }
        }
        // Expand the attraction regions of local minimizers until they do not overlap
        for i in 0..num_minima {
            if i != 1 {
                temp_min = f64::MAX;
                for j in 0..num_minima {
                    if i != j {
                        let temp_d1 =
                            norm(&minima.local_min[i], &minima.local_min[j]) - minima.rho[j];
                        if temp_d1 < temp_min {
                            temp_min = temp_d1;
                        }
                    }
                }
                if temp_min > minima.rho[i] + options.precision {
                    minima.rho[i] = temp_min;
                }
            }
        }
        // Correct the radii by weight coefficients w(i)
        for i in 0..num_minima {
            minima.rho[i] *= minima.w_rho[i];
        }
        // Set the local minima values f(i) of test functions
        minima.peak[0] = 0.0;
        minima.peak[1] = 0.0;
        for i in 2..num_minima {
            let random = rng.gen();
            temp_d1 = norm(&minima.local_min[0], &minima.local_min[i]);
            temp_min = (minima.rho[i] - temp_d1).mul_add(minima.rho[i] - temp_d1, minima.f[0]);
            temp_d1 = (1.0 + random) * minima.rho[i];
            temp_d2 = random * (temp_min - global_value);
            if temp_d2 < temp_d1 {
                temp_d1 = temp_d2;
            }
            minima.peak[i] = temp_d1;
            minima.f[i] = temp_min - minima.peak[i];
        }
        // Find all possible global minimizers and create a list of their indices among all the minimizers
        glob.num_global_minima = 0;
        for i in 0..num_minima {
            if (minima.f[i] >= global_value - options.precision)
                && (minima.f[i] <= global_value + options.precision)
            {
                glob.gm_index[glob.num_global_minima] = i;
                glob.num_global_minima += 1;
            } else {
                glob.gm_index[num_minima - 1 - i + glob.num_global_minima] = i;
            }
        }
        if glob.num_global_minima == 0 {
            return Err(GKLSError::FloatingPoint);
        }
        Ok(())
    }

    /// Generate a new test function
    pub fn new(
        nf: usize,
        options: Options,
        dim: usize,
        num_minima: usize,
        global_value: f64,
        global_radius: f64,
        global_dist: f64,
    ) -> Result<Self, GKLSError> {
        let mut sin_phi: f64;
        let gap = global_radius;
        if !(1..=100).contains(&nf) {
            return Err(GKLSError::FuncNumber);
        }
        let domain = Domain::default(dim);
        Self::parameters_check(
            dim,
            num_minima,
            &domain,
            &options,
            global_value,
            global_radius,
            global_dist,
        )?;
        let mut minima: Minima = Minima {
            local_min: vec![vec![0.0; dim]; num_minima],
            f: vec![0.0; num_minima],
            w_rho: vec![0.0; num_minima],
            peak: vec![0.0; num_minima],
            rho: vec![0.0; num_minima],
        };
        let seed = (nf - 1) + (num_minima - 1) * 100 + dim * 1_000_000;
        let mut rng = ranf::Ranf::new(100, 37, 70, 1009, 1009, seed.try_into().unwrap());
        for i in 0..dim {
            let r = rng.gen();
            minima.local_min[0][i] = r.mul_add(domain.right[i] - domain.left[i], domain.left[i]);
        }
        rng.reset();
        minima.f[0] = options.paraboloid_min;
        let w_phi = PI * rng.gen();
        minima.local_min[1][0] = global_dist.mul_add(f64::cos(w_phi), minima.local_min[0][0]);
        if minima.local_min[1][0] > domain.right[0] - options.precision
            || minima.local_min[1][0] < domain.left[0] + options.precision
        {
            minima.local_min[1][0] = global_dist.mul_add(-f64::cos(w_phi), minima.local_min[0][0]);
        }
        sin_phi = f64::sin(w_phi);
        for j in 1..dim - 1 {
            let w_2 = PI * rng.gen();
            minima.local_min[1][j] =
                (global_dist * f64::cos(2.0 * w_2)).mul_add(sin_phi, minima.local_min[0][j]);
            if minima.local_min[1][j] > domain.right[j] - options.precision
                || minima.local_min[1][j] < domain.left[j] + options.precision
            {
                minima.local_min[1][j] =
                    (global_dist * f64::cos(2.0 * w_2)).mul_add(-sin_phi, minima.local_min[0][j]);
            }
            sin_phi *= f64::sin(2.0 * w_2);
        }
        minima.local_min[1][dim - 1] = global_dist.mul_add(sin_phi, minima.local_min[0][dim - 1]);
        if minima.local_min[1][dim - 1] > domain.right[dim - 1] - options.precision
            || minima.local_min[1][dim - 1] < domain.left[dim - 1] + options.precision
        {
            minima.local_min[1][dim - 1] =
                global_dist.mul_add(-sin_phi, minima.local_min[0][dim - 1]);
        }
        minima.f[1] = global_value;
        for i in 0..num_minima {
            minima.w_rho[i] = 0.99;
        }
        minima.w_rho[1] = 1.0;
        let delta = options.delta_max_value * rng.gen();
        while let CoincidenceCondition::LocalMinCoincidence
        | CoincidenceCondition::ParabolaMinCoincidence =
            coincidence_check(&minima, num_minima, options.precision)
        {
            let mut i = 2;
            while i < num_minima {
                loop {
                    rng.reset();
                    for j in 0..dim {
                        minima.local_min[i][j] = rng
                            .gen()
                            .mul_add(domain.right[j] - domain.left[j], domain.left[j]);
                    }
                    if (global_radius + gap) - norm(&minima.local_min[i], &minima.local_min[1])
                        < options.precision
                    {
                        break;
                    }
                }
                i += 1;
            }
        }
        let mut glob = GlobalMinima {
            num_global_minima: 0,
            gm_index: vec![0; num_minima],
        };
        Self::set_basins(
            num_minima,
            &mut minima,
            &mut glob,
            global_value,
            global_radius,
            &options,
            rng,
        )?;
        Ok(Self {
            minima,
            domain,
            num_minima,
            _global_value: global_value,
            _global_radius: global_radius,
            _global_dist: global_dist,
            delta,
            options,
            dim,
            _glob: glob,
        })
    }

    /// Evaluates the generated function of the ND-type (non-differentiable)
    #[must_use]
    pub fn nd_func(&self, x: &[f64]) -> f64 {
        let mut index: usize = 0;
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        while index < self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        let mut norm_: f64;
        if index == self.num_minima {
            norm_ = norm(&self.minima.local_min[0], x);
            return norm_.mul_add(norm_, self.minima.f[0]);
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return self.minima.f[index];
        }
        norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_.mul_add(norm_, self.minima.f[0]) - self.minima.f[index];
        let rho = self.minima.rho[index];
        norm_ = norm(&self.minima.local_min[index], x);
        let mut scal = 0.0;
        for (val, local_min0, local_min_index) in izip!(
            x,
            self.minima.local_min[0].iter(),
            self.minima.local_min[index].iter()
        ) {
            scal += (val - local_min_index) * (local_min0 - local_min_index);
        }
        (1.0 - 2.0 / rho * scal / norm_ + a / rho.powi(2))
            .mul_add(norm_.powi(2), self.minima.f[index])
    }

    /// Evaluates the generated function of the D-type (continuously differentiable)
    #[must_use]
    pub fn d_func(&self, x: &[f64]) -> f64 {
        let mut norm_: f64;
        let mut scal: f64;

        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index: usize = 1;
        while index < self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        if index == self.num_minima {
            norm_ = norm(&self.minima.local_min[0], x);
            return norm_.mul_add(norm_, self.minima.f[0]);
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return self.minima.f[index];
        }
        norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_.mul_add(norm_, self.minima.f[0]) - self.minima.f[index];
        let rho = self.minima.rho[index];
        norm_ = norm(&self.minima.local_min[index], x);
        scal = 0.0;
        for (val, local_min0, local_min_index) in izip!(
            x,
            self.minima.local_min[0].iter(),
            self.minima.local_min[index].iter()
        ) {
            scal += (val - local_min_index) * (local_min0 - local_min_index);
        }
        (2.0 / rho.powi(2) * scal / norm_ - 2.0 * a / rho.powi(3)).mul_add(
            norm_.powi(3),
            (1.0 - 4.0 * scal / norm_ / rho + 3.0 * a / rho.powi(2)) * norm_.powi(2),
        ) + self.minima.f[index]
    }

    /// Evaluates the generated function of the D2-type (twice continuously differentiable)
    #[must_use]
    pub fn d2_func(&self, x: &[f64]) -> f64 {
        let mut norm_: f64;
        let mut scal: f64;
        let mut index: usize = 1;
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        while index < self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        if index == self.num_minima {
            norm_ = norm(&self.minima.local_min[0], x);
            return norm_.mul_add(norm_, self.minima.f[0]);
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return self.minima.f[index];
        }
        norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_.mul_add(norm_, self.minima.f[0]) - self.minima.f[index];
        let rho = self.minima.rho[index];
        norm_ = norm(&self.minima.local_min[index], x);
        scal = 0.0;
        for (val, local_min0, local_min_index) in izip!(
            x,
            self.minima.local_min[0].iter(),
            self.minima.local_min[index].iter()
        ) {
            scal += (val - local_min_index) * (local_min0 - local_min_index);
        }
        let term1: f64 = (-6.0 * scal / norm_ / rho + 6.0 * a / rho.powi(2) + 1.0
            - self.delta / 2.0)
            * norm_.powi(2)
            / rho.powi(2);
        let term2: f64 = 1.5f64.mul_add(
            self.delta,
            16.0 * scal / norm_ / rho - 15.0 * a / rho / rho - 3.0,
        ) * norm_
            / rho;
        let term3: f64 = 1.5f64.mul_add(
            -self.delta,
            -12.0 * scal / norm_ / rho + 10.0 * a / rho.powi(2) + 3.0,
        );
        let term4: f64 = 0.5 * self.delta * norm_ * norm_;
        (term1 + term2 + term3) * norm_.powi(3) / rho + term4 + self.minima.f[index]
    }

    /// Evaluates the first order partial derivative of D-type function with respet to the
    /// variable wiht index `var_j`
    #[must_use]
    pub fn d_deriv(&self, var_j: usize, x: &[f64]) -> f64 {
        let mut var_j = var_j;
        if var_j == 0 || var_j > self.dim {
            return self.options.max_value;
        }
        var_j -= 1;
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index = 1;
        while index < self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        if index == self.num_minima {
            return 2.0 * (x[var_j] - self.minima.local_min[0][var_j]);
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return 0.0;
        }
        let mut norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_.mul_add(norm_, self.minima.f[0]) - self.minima.f[index];
        let rho = self.minima.rho[index];
        norm_ = norm(&self.minima.local_min[index], x);
        let mut scal = 0.0;
        for (val, local_min0, local_min_index) in izip!(
            x,
            self.minima.local_min[0].iter(),
            self.minima.local_min[index].iter()
        ) {
            scal += (*val - local_min_index) * (local_min0 - local_min_index);
        }
        let dif = x[var_j] - self.minima.local_min[index][var_j];
        let h = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j])
            .mul_add(norm_, -scal * dif / norm_);
        h.mul_add(
            (2.0 / rho.powi(2)).mul_add(norm_, -4.0 / rho),
            dif * ((6.0 / rho.powi(2)).mul_add(
                a,
                (8.0 / rho / norm_).mul_add(
                    -scal,
                    (6.0 / rho.powi(2)).mul_add(scal, -6.0 / rho.powi(3) * a * norm_),
                ),
            ) + 2.0),
        )
    }

    /// Evaluates the first order partial derivative of D2-type function with respect
    /// to the variable with index `var_j`
    #[must_use]
    pub fn d2_deriv1(&self, var_j: usize, x: &[f64]) -> f64 {
        let mut var_j = var_j;
        if var_j == 0 || var_j > self.dim {
            return self.options.max_value;
        }
        var_j -= 1;
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index = 1;
        while index < self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        if index == self.num_minima {
            return 2.0 * (x[var_j] - self.minima.local_min[0][var_j]);
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return 0.0;
        }
        let mut norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_.mul_add(norm_, self.minima.f[0]) - self.minima.f[index];
        let rho = self.minima.rho[index];
        norm_ = norm(&self.minima.local_min[index], x);
        let mut scal = 0.0;
        for (val, local_min0, local_min_index) in izip!(
            x,
            self.minima.local_min[0].iter(),
            self.minima.local_min[index].iter()
        ) {
            scal += (*val - local_min_index) * (local_min0 - local_min_index);
        }
        let dif = x[var_j] - self.minima.local_min[index][var_j];
        let h = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j])
            .mul_add(norm_, -scal * dif / norm_);
        dif.mul_add(
            self.delta,
            (h * norm_ / rho.powi(2)).mul_add(
                -6.0 * norm_.powi(2) / rho.powi(2) + 16.0 * norm_ / rho - 12.0,
                dif * norm_
                    * ((2.5f64.mul_add(
                        -self.delta,
                        (-30.0 / rho / norm_).mul_add(scal, 30.0 / rho.powi(2) * a) + 5.0,
                    ) / rho.powi(3))
                    .mul_add(
                        norm_.powi(2),
                        6.0f64.mul_add(
                            self.delta,
                            (64.0 / rho / norm_).mul_add(scal, -60.0 / rho.powi(2) * a) - 12.0,
                        ) / rho.powi(2)
                            * norm_,
                    ) + 4.5f64.mul_add(
                        -self.delta,
                        (-36.0 / rho / norm_).mul_add(scal, 30.0 / rho.powi(2) * a) + 9.0,
                    ) / rho),
            ),
        )
    }

    /// Evaluates the second order partial derivative of D2-type function with respect
    /// to the variables with indexes `var_j` and `var_k`
    #[must_use]
    pub fn d2_deriv2(&self, var_j: usize, var_k: usize, x: &[f64]) -> f64 {
        let mut var_j = var_j;
        let mut var_k = var_k;
        if var_j == 0 || var_j > self.dim {
            return self.options.max_value;
        }
        if var_k == 0 || var_k > self.dim {
            return self.options.max_value;
        }
        let the_same = var_j == var_k;
        var_j -= 1;
        var_k -= 1;
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index = 1;
        while index < self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        if index == self.num_minima {
            if the_same {
                return 2.0;
            }
            return 0.0;
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            if the_same {
                return self.delta;
            }
            return 0.0;
        }
        let mut norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_.mul_add(norm_, self.minima.f[0]) - self.minima.f[index];
        let rho = self.minima.rho[index];
        norm_ = norm(&self.minima.local_min[index], x);
        let mut scal = 0.0;
        for (val, local_min0, local_min_index) in izip!(
            x,
            self.minima.local_min[0].iter(),
            self.minima.local_min[index].iter()
        ) {
            scal += (val - local_min_index) * (local_min0 - local_min_index);
        }
        let difj = x[var_j] - self.minima.local_min[index][var_j];
        let difk = x[var_k] - self.minima.local_min[index][var_k];
        let hj = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j])
            .mul_add(norm_, -scal * difj / norm_);
        let hk = (self.minima.local_min[0][var_k] - self.minima.local_min[index][var_k])
            .mul_add(norm_, -scal * difk / norm_);
        let mut dh = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j]) * difk
            / norm_
            - hk * difj / norm_.powi(2);
        if the_same {
            dh -= scal / norm_;
        }
        let mut dq_jk = (12.0 / rho.powi(2)).mul_add(
            -dh.mul_add(norm_, hj * difk / norm_),
            (8.0 / rho.powi(2)
                * 1.5f64.mul_add(
                    self.delta,
                    (15.0 / rho.powi(2)).mul_add(-a, 16.0 / rho * scal / norm_) - 3.0,
                )
                * difj)
                .mul_add(
                    difk,
                    (64.0 / rho.powi(3) * hk).mul_add(
                        difj,
                        (16.0 / rho.powi(3)).mul_add(
                            dh.mul_add(norm_.powi(2), 2.0 * hj * difk),
                            (15.0 / rho.powi(3)
                                * 0.5f64.mul_add(
                                    -self.delta,
                                    (6.0 / rho.powi(2)).mul_add(a, -6.0 / rho * scal / norm_) + 1.0,
                                )
                                * difj
                                * difk)
                                .mul_add(
                                    norm_,
                                    (-6.0 / rho.powi(4)).mul_add(
                                        dh.mul_add(norm_.powi(3), 3.0 * hj * difk * norm_),
                                        -30.0 / rho.powi(4) * hk * difj * norm_,
                                    ),
                                ),
                        ),
                    ),
                ),
        ) - 36.0 / rho.powi(2) * hk * difj / norm_
            + 3.0 / rho
                * 1.5f64.mul_add(
                    -self.delta,
                    (10.0 / rho.powi(2)).mul_add(a, -12.0 / rho * scal / norm_) + 3.0,
                )
                * difj
                * difk
                / norm_;
        if the_same {
            dq_jk = (norm_ / rho).mul_add(
                1.5f64.mul_add(
                    -self.delta,
                    (10.0 / rho.powi(2)).mul_add(a, -12.0 / rho * scal / norm_) + 3.0,
                ),
                (4.0 * norm_.powi(2) / rho.powi(2)).mul_add(
                    1.5f64.mul_add(
                        self.delta,
                        (15.0 / rho.powi(2)).mul_add(-a, 16.0 / rho * scal / norm_) - 3.0,
                    ),
                    (5.0 * norm_.powi(3) / rho.powi(3)).mul_add(
                        0.5f64.mul_add(
                            -self.delta,
                            (6.0 / rho.powi(2)).mul_add(a, -6.0 / rho * scal / norm_) + 1.0,
                        ),
                        dq_jk,
                    ),
                ),
            ) + self.delta;
        }
        dq_jk
    }

    /// Evaluates the gradient of the D-type function
    #[must_use]
    pub fn d_gradient(&self, x: &[f64]) -> Option<Vec<f64>> {
        let mut g = vec![f64::NAN; self.dim];
        for (i, g_) in g.iter_mut().enumerate() {
            *g_ = self.d_deriv(i + 1, x);
            if *g_ >= self.options.max_value - 1000.0 {
                return None;
            }
        }
        Some(g)
    }

    /// Evaluates the gradient of the D2-type function
    #[must_use]
    pub fn d2_gradient(&self, x: &[f64]) -> Option<Vec<f64>> {
        let mut g = vec![f64::NAN; self.dim];
        for (i, g_) in g.iter_mut().enumerate() {
            *g_ = self.d2_deriv1(i + 1, x);
            if *g_ >= self.options.max_value - 1000.0 {
                return None;
            }
        }
        Some(g)
    }

    /// Evaluates the hessian of the D2-type function
    #[must_use]
    pub fn d2_hessian(&self, x: &[f64]) -> Option<Vec<Vec<f64>>> {
        let mut h = vec![vec![f64::NAN; self.dim]; self.dim];
        for (i, h_) in h.iter_mut().enumerate() {
            for (j, h__) in h_.iter_mut().enumerate() {
                *h__ = self.d2_deriv2(i + 1, j + 1, x);
                if *h__ >= self.options.max_value - 1000.0 {
                    return None;
                }
            }
        }
        Some(h)
    }
}

#[cfg(feature = "test_cbinding")]
#[cfg(test)]
mod tests {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    const SEED: u64 = 42;
    const N_TESTPOINTS: usize = 10;
    const TOLERANCE: f64 = 10e-6;

    fn generate_testpoints(dim: usize) -> Vec<Vec<f64>> {
        let mut rng = StdRng::seed_from_u64(SEED);
        (0..N_TESTPOINTS)
            .map(|_| (0..dim).map(|_| rng.gen()).collect())
            .collect()
    }

    macro_rules! test_problem_value_function {
        ($gkls_problem:expr, $gkls_problem_cimpl:expr, $method:ident, $tolerance:expr) => {
            #[test]
            fn $method() {
                use crate::tests::generate_testpoints;
                let tolerance = $tolerance;
                let problem = $gkls_problem;
                let problem_cimpl = $gkls_problem_cimpl;
                for input in generate_testpoints(problem.dim) {
                    {
                        let result = problem.$method(&input);
                        let result_c_impl = problem_cimpl.$method(&input);
                        assert!(
                            (result - result_c_impl).abs() <= tolerance,
                            "{} is not approximately equal to \n{} with absolute tolerance {}",
                            result,
                            result_c_impl,
                            tolerance
                        );
                        if result_c_impl != 0.0 {
                            assert!(
                                ((result - result_c_impl).abs()) / result_c_impl <= tolerance,
                                "{} is not approximately equal to \n{} with relative tolerance {}",
                                result,
                                result_c_impl,
                                tolerance
                            );
                        }
                    }
                }
            }
        };
    }

    macro_rules! test_problem_deriv1_function {
        ($gkls_problem:expr, $gkls_problem_cimpl:expr, $method:ident, $tolerance:expr) => {
            #[test]
            fn $method() {
                use crate::tests::generate_testpoints;
                let tolerance = $tolerance;
                let problem = $gkls_problem;
                let problem_cimpl = $gkls_problem_cimpl;
                for input in generate_testpoints(problem.dim) {
                    for dimension in 0..problem.dim {
                        let result = problem.$method(dimension, &input);
                        let result_c_impl = problem_cimpl.$method(dimension, &input);
                        assert!(
                            (result - result_c_impl).abs() <= tolerance,
                            "{} is not approximately equal to \n{} with absolute tolerance {}",
                            result,
                            result_c_impl,
                            tolerance
                        );
                        if result_c_impl != 0.0 {
                            assert!(
                                ((result - result_c_impl).abs()) / result_c_impl <= tolerance,
                                "{} is not approximately equal to \n{} with relative tolerance {}",
                                result,
                                result_c_impl,
                                tolerance
                            );
                        }
                    }
                }
            }
        };
    }

    macro_rules! test_problem_deriv2_function {
        ($gkls_problem:expr, $gkls_problem_cimpl:expr, $tolerance:expr) => {
            #[test]
            fn d2_deriv2() {
                use crate::tests::generate_testpoints;
                let tolerance = $tolerance;
                let problem = $gkls_problem;
                let problem_cimpl = $gkls_problem_cimpl;
                for input in generate_testpoints(problem.dim) {
                    for dim0 in 0..problem.dim {
                        for dim1 in 0..problem.dim {
                            let result = problem.d2_deriv2(dim0, dim1, &input);
                            let result_c_impl = problem_cimpl.d2_deriv2(dim0, dim1, &input);
                            assert!(
                                (result - result_c_impl).abs() <= tolerance,
                                "{} is not approximately equal to \n{} with absolute tolerance {}",
                                result,
                                result_c_impl,
                                tolerance
                            );
                        if result_c_impl != 0.0 {
                            assert!(
                                ((result - result_c_impl).abs()) / result_c_impl <= tolerance,
                                "{} is not approximately equal to \n{} with relative tolerance {}",
                                result,
                                result_c_impl,
                                tolerance
                            );
                        }
                        }
                    }
                }
            }
        };
    }
    macro_rules! test_problem_gradient {
        ($gkls_problem:expr, $gkls_problem_cimpl:expr, $method:ident, $tolerance:expr) => {
            #[test]
            fn $method() {
                use crate::tests::generate_testpoints;
                let tolerance = $tolerance;
                let problem = $gkls_problem;
                let problem_cimpl = $gkls_problem_cimpl;
                for input in generate_testpoints(problem.dim) {
                    let result = problem.$method(&input);
                    let result_c_impl = problem_cimpl.$method(&input);
                    match (result, result_c_impl) {
                        (Some(result), Some(result_c_impl)) => {
                        for i in 0..result.len() {
                            assert!(
                                (result[i] - result_c_impl[i]).abs() <= tolerance,
                                "{} is not approximately equal to \n{} with absolute tolerance {}",
                                result[i],
                                result_c_impl[i],
                                tolerance
                            );
                            if result_c_impl[i] != 0.0 {
                                assert!(
                                    ((result[i] - result_c_impl[i]).abs()) / result_c_impl[i]
                                        <= tolerance,
                                    "{} is not approximately equal to \n{} with relative tolerance {}",
                                    result[i],
                                    result_c_impl[i],
                                    tolerance
                                );
                            }
                            }
                        },
                        (Some(_), None) => panic!("Rust implementation returned values while C implementation returned None"),
                        (None, Some(_)) => panic!("C implementation returned values while Rust implementation returned None"),
                        (None, None) => ()
                        // _ => ()
                    }
                }
            }
        };
    }
    macro_rules! test_problem_hessian {
        ($gkls_problem:expr, $gkls_problem_cimpl:expr, $tolerance:expr) => {
            #[test]
            fn hessian() {
                use crate::tests::generate_testpoints;
                let tolerance = $tolerance;
                let problem = $gkls_problem;
                let problem_cimpl = $gkls_problem_cimpl;
                for input in generate_testpoints(problem.dim) {
                    let result = problem.d2_hessian(&input);
                    let result_c_impl = problem_cimpl.d2_hessian(&input);
                    match (result, result_c_impl) {
                        (Some(result), Some(result_c_impl)) => {
                        for i in 0..result.len() {
                            for j in 0..result[i].len() {
                                assert!(
                                    (result[i][j] - result_c_impl[i][j]).abs() <= tolerance,
                                    "{} is not approximately equal to \n{} with absolute tolerance {}",
                                    result[i][j],
                                    result_c_impl[i][j],
                                    tolerance
                                );
                            if result_c_impl[i][j] != 0.0 {
                                assert!(
                                    ((result[i][j] - result_c_impl[i][j]).abs()) / result_c_impl[i][j]
                                        <= tolerance,
                                    "{} is not approximately equal to \n{} with relative tolerance {}",
                                    result[i][j],
                                    result_c_impl[i][j],
                                    tolerance
                                );
                            }
                            }
                        }
                        },
                        (Some(_), None) => panic!("Rust implementation returned values while C implementation returned None"),
                        (None, Some(_)) => panic!("C implementation returned values while Rust implementation returned None"),
                        (None, None) => ()
                        // _ => ()
                    }
                }
            }
        };
    }

    mod default_problem {
        use crate::{c_gkls::CGKLSProblem, tests::TOLERANCE, Problem};
        test_problem_value_function!(
            Problem::default(),
            CGKLSProblem::default(),
            nd_func,
            TOLERANCE
        );
        test_problem_value_function!(
            Problem::default(),
            CGKLSProblem::default(),
            d_func,
            TOLERANCE
        );
        test_problem_value_function!(
            Problem::default(),
            CGKLSProblem::default(),
            d2_func,
            TOLERANCE
        );
        test_problem_deriv1_function!(
            Problem::default(),
            CGKLSProblem::default(),
            d_deriv,
            TOLERANCE
        );
        test_problem_deriv1_function!(
            Problem::default(),
            CGKLSProblem::default(),
            d2_deriv1,
            TOLERANCE
        );
        test_problem_deriv2_function!(Problem::default(), CGKLSProblem::default(), TOLERANCE);
    }
    macro_rules! test_set {
        ($name:ident, $nf:expr, $dim:expr, $num_minima:expr, $global_value:expr, $global_radius:expr, $global_dist:expr) => {
            mod $name {
                use crate::tests::TOLERANCE;
                use crate::{c_gkls::CGKLSProblem, Options, Problem};
                test_problem_value_function!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    nd_func,
                    TOLERANCE
                );
                test_problem_value_function!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    d_func,
                    TOLERANCE
                );
                test_problem_value_function!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    d2_func,
                    TOLERANCE
                );
                test_problem_deriv1_function!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    d_deriv,
                    TOLERANCE
                );
                test_problem_deriv1_function!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    d2_deriv1,
                    TOLERANCE
                );
                test_problem_deriv2_function!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    TOLERANCE
                );
                test_problem_gradient!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    d_gradient,
                    TOLERANCE
                );
                test_problem_gradient!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    d2_gradient,
                    TOLERANCE
                );
                test_problem_hessian!(
                    Problem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    CGKLSProblem::new(
                        $nf,
                        Options::default(),
                        $dim,
                        $num_minima,
                        $global_value,
                        $global_radius,
                        $global_dist
                    )
                    .unwrap(),
                    TOLERANCE
                );
            }
        };
    }

    test_set!(nf3dim7, 3, 7, 5, -1.0, 1. / 3., 2. / 3.);
    test_set!(nf7dim7, 7, 7, 5, -1.0, 1. / 3., 2. / 3.);
    test_set!(nf42dim13, 42, 13, 5, -1.0, 1. / 3., 2. / 3.);
    test_set!(nf17dim71, 17, 71, 5, -1.7, 2. / 9., 7. / 8.);
}
