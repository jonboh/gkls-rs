use std::f64::consts::PI;

use itertools::izip;
use rand::Rng;

pub struct Options {
    max_value: f64,
    precision: f64,
    paraboloid_min: f64,
    global_min_value: f64,
    delta_max_value: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_value: 1e+100,
            precision: 1e-10,
            paraboloid_min: 0.0,
            global_min_value: -1.0,
            delta_max_value: 10.0,
        }
    }
}

#[derive(Debug)]
pub enum GKLSError {
    Dim,
    NumMinima,
    FuncNumber,
    Boundary,
    GlobalMinValue,
    GlobalDist,
    GlobalRadius,
    Memory,
    DerivEval,
    FloatingPoint,
}

struct Minima {
    local_min: Vec<Vec<f64>>,
    f: Vec<f64>,
    w_rho: Vec<f64>,
    peak: Vec<f64>,
    rho: Vec<f64>,
}

struct GlobalMinima {
    num_global_minima: usize,
    gm_index: Vec<usize>,
}

enum CoincidenceCondition {
    ParabolaMinCoincidence,
    LocalMinCoincidence,
    Ok,
}

pub struct Domain {
    left: Vec<f64>,
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

pub struct Problem {
    options: Options,
    domain: Domain,
    dim: usize,
    num_minima: usize,
    global_dist: f64,
    global_radius: f64,
    global_value: f64,
    delta: f64,
    minima: Minima,
    glob: GlobalMinima,
}

impl Default for Problem {
    fn default() -> Self {
        let nf = 0;
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
    // TODO: probably don't need num_minima and can len() minima
    let num_minima = num_minima;
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
        // NOTE: maybe not pub, or removed
        if dim <= 1 {
            return Err(GKLSError::Dim);
        }
        if num_minima <= 1 {
            return Err(GKLSError::NumMinima);
        }
        let dim = dim;
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
        if global_dist >= 0.5 * min_side - options.precision || global_dist <= options.precision {
            return Err(GKLSError::GlobalDist);
        }
        if global_radius >= 0.5 * global_dist + options.precision
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
        global_dist: f64,
        global_radius: f64,
        options: &Options,
    ) -> Result<(), GKLSError> {
        let mut temp_min: f64;
        let mut temp_d1: f64;
        let mut temp_d2: f64;
        let mut dist: f64;
        let num_minima = num_minima;
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
            minima.rho[i] = minima.w_rho[i] * minima.rho[i];
        }
        // Set the local minima values f(i) of test functions
        minima.peak[0] = 0.0;
        minima.peak[1] = 0.0;
        let mut rng = rand::thread_rng();
        for i in 2..num_minima {
            temp_d1 = norm(&minima.local_min[0], &minima.local_min[i]);
            temp_min = (minima.rho[i] - temp_d1).powi(2) + minima.f[0];
            temp_d1 = (1.0 + rng.gen::<f64>()) * minima.rho[i];
            temp_d2 = rng.gen::<f64>() * (temp_min - global_value);
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

    pub fn new(
        nf: usize,
        options: Options,
        dim: usize,
        num_minima: usize,
        global_value: f64,
        global_radius: f64,
        global_dist: f64,
    ) -> Result<Problem, GKLSError> {
        let mut error: i32;
        let mut i: usize;
        let mut j: usize;
        let mut sin_phi: f64;
        let mut gap = global_radius;
        if nf < 1 || nf > 100 {
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
        let mut rng = rand::thread_rng();
        for i in 0..dim {
            minima.local_min[0][i] =
                domain.left[i] + rng.gen::<f64>() * (domain.right[i] - domain.left[i]);
        }
        minima.f[0] = options.paraboloid_min;
        minima.local_min[1][0] =
            minima.local_min[0][0] + global_dist * f64::cos(PI * rng.gen::<f64>());
        if minima.local_min[1][0] > domain.right[0] - options.precision
            || minima.local_min[1][0] < domain.left[0] + options.precision
        {
            minima.local_min[1][0] =
                minima.local_min[0][0] - global_dist * f64::cos(PI * rng.gen::<f64>());
        }
        sin_phi = f64::sin(PI * rng.gen::<f64>());
        for j in 1..dim - 1 {
            minima.local_min[1][j] = minima.local_min[0][j]
                + global_dist * f64::cos(2.0 * PI * rng.gen::<f64>()) * sin_phi;
            if minima.local_min[1][j] > domain.right[j] - options.precision
                || minima.local_min[1][j] < domain.left[j] + options.precision
            {
                minima.local_min[1][j] = minima.local_min[0][j]
                    - global_dist * f64::cos(2.0 * PI * rng.gen::<f64>()) * sin_phi;
            }
            sin_phi *= f64::sin(2.0 * PI * rng.gen::<f64>());
        }
        minima.local_min[1][dim - 1] = minima.local_min[0][dim - 1] + global_dist * sin_phi;
        if minima.local_min[1][dim - 1] > domain.right[dim - 1] - options.precision
            || minima.local_min[1][dim - 1] < domain.left[dim - 1] + options.precision
        {
            minima.local_min[1][dim - 1] = minima.local_min[0][dim - 1] - global_dist * sin_phi;
        }
        minima.f[1] = global_value;
        for i in 0..num_minima {
            minima.w_rho[i] = 0.99;
        }
        minima.w_rho[1] = 1.0;
        let delta = options.delta_max_value * rng.gen::<f64>();
        let mut i = 2;
        while let CoincidenceCondition::LocalMinCoincidence
        | CoincidenceCondition::ParabolaMinCoincidence =
            coincidence_check(&minima, num_minima, options.precision)
        {
            while i < num_minima {
                loop {
                    for j in 0..dim {
                        minima.local_min[i][j] =
                            domain.left[j] + rng.gen::<f64>() * (domain.right[j] - domain.left[j]);
                    }
                    if global_radius + gap - norm(&minima.local_min[i], &minima.local_min[1])
                        > options.precision
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
            global_dist,
            global_radius,
            &options,
        )?;
        Ok(Self {
            minima,
            domain,
            num_minima,
            global_value,
            global_radius,
            global_dist,
            delta,
            options,
            dim,
            glob,
        })
    }

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
            return norm_ * norm_ + self.minima.f[0];
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return self.minima.f[index];
        }
        norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_ * norm_ + self.minima.f[0] - self.minima.f[index];
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
        (1.0 - 2.0 / rho * scal / norm_ + a / rho / rho) * norm_ * norm_ + self.minima.f[index]
    }

    pub fn d_func(&self, x: &[f64]) -> f64 {
        let mut norm_: f64;
        let mut scal: f64;

        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index: usize = 0;
        while index <= self.num_minima
            && norm(&self.minima.local_min[index], x) > self.minima.rho[index]
        {
            index += 1;
        }
        if index == self.num_minima {
            norm_ = norm(&self.minima.local_min[0], x);
            return norm_ * norm_ + self.minima.f[0];
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return self.minima.f[index];
        }
        norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_ * norm_ + self.minima.f[0] - self.minima.f[index];
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
        (2.0 / rho / rho * scal / norm_ - 2.0 * a / rho / rho / rho) * norm_.powi(3)
            + (1.0 - 4.0 * scal / norm_ / rho + 3.0 * a / rho / rho) * norm_.powi(2)
            + self.minima.f[index]
    }

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
            return norm_ * norm_ + self.minima.f[0];
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            return self.minima.f[index];
        }
        norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_ * norm_ + self.minima.f[0] - self.minima.f[index];
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
        let delta: f64 = 1.0 - 1.0 / (2.0 * self.dim as f64);
        let term1: f64 = (-6.0 * scal / norm_ / rho + 6.0 * a / rho / rho + 1.0 - delta / 2.0)
            * norm_.powi(2)
            / rho.powi(2);
        let term2: f64 =
            (16.0 * scal / norm_ / rho - 15.0 * a / rho / rho - 3.0 + 1.5 * delta) * norm_ / rho;
        let term3: f64 = (-12.0 * scal / norm_ / rho + 10.0 * a / rho / rho + 3.0 - 1.5 * delta)
            * norm_.powi(3)
            / rho;
        let term4: f64 = 0.5 * delta * norm_.powi(2);
        term1 + term2 + term3 + term4 + self.minima.f[index]
    }

    fn d_deriv(&self, var_j: usize, x: &[f64]) -> f64 {
        if var_j == 0 || var_j >= self.dim {
            return self.options.max_value;
        }
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index = 0;
        while index <= self.num_minima
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
        let a = norm_ * norm_ + self.minima.f[0] - self.minima.f[index];
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
        let h = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j]) * norm_
            - scal * dif / norm_;
        h * (2.0 / rho / rho * norm_ - 4.0 / rho)
            + dif
                * (6.0 / rho / rho * scal
                    - 6.0 / rho / rho / rho * a * norm_
                    - 8.0 / rho / norm_ * scal
                    + 6.0 / rho / rho * a
                    + 2.0)
    }

    fn d2_deriv1(&self, var_j: usize, x: &[f64]) -> f64 {
        if var_j == 0 || var_j >= self.dim {
            return self.options.max_value;
        }
        for (val, left, right) in izip!(x, self.domain.left.iter(), self.domain.right.iter()) {
            if *val < left - self.options.precision || *val > right + self.options.precision {
                return self.options.max_value;
            }
        }
        let mut index = 0;
        while index <= self.num_minima
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
        let a = norm_ * norm_ + self.minima.f[0] - self.minima.f[index];
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
        let h = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j]) * norm_
            - scal * dif / norm_;
        let delta: f64 = 1.0;
        h * norm_ / rho / rho * (-6.0 * norm_ * norm_ / rho / rho + 16.0 * norm_ / rho - 12.0)
            + dif
                * norm_
                * ((-30.0 / rho / norm_ * scal + 30.0 / rho / rho * a + 5.0 - 2.5 * delta)
                    / rho
                    / rho
                    / rho
                    * norm_
                    * norm_
                    + (64.0 / rho / norm_ * scal - 60.0 / rho / rho * a - 12.0 + 6.0 * delta)
                        / rho
                        / rho
                        * norm_
                    + (-36.0 / rho / norm_ * scal + 30.0 / rho / rho * a + 9.0 - 4.5 * delta) / rho)
            + dif * delta
    }

    fn d2_deriv2(&self, var_j: usize, var_k: usize, x: &[f64]) -> f64 {
        if var_j == 0 || var_j >= self.dim {
            return self.options.max_value;
        }
        if var_k == 0 || var_k >= self.dim {
            return self.options.max_value;
        }
        let the_same = var_j == var_k;
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
            } else {
                return 0.0;
            }
        }
        if norm(x, &self.minima.local_min[index]) < self.options.precision {
            if the_same {
                return self.delta;
            } else {
                return 0.0;
            }
        }
        let mut norm_ = norm(&self.minima.local_min[0], &self.minima.local_min[index]);
        let a = norm_ * norm_ + self.minima.f[0] - self.minima.f[index];
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
        let hj = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j]) * norm_
            - scal * difj / norm_;
        let hk = (self.minima.local_min[0][var_k] - self.minima.local_min[index][var_k]) * norm_
            - scal * difk / norm_;
        let mut dh = (self.minima.local_min[0][var_j] - self.minima.local_min[index][var_j]) * difk
            / norm_
            - hk * difj / norm_.powi(2);
        if the_same {
            dh -= scal / norm_;
        }
        let mut dQ_jk = -6.0 / rho.powi(4) * (dh * norm_.powi(3) + 3.0 * hj * difk * norm_.powi(2))
            - 30.0 / rho.powi(4) * hk * difj * norm_
            + 15.0 / rho.powi(3)
                * (-6.0 / rho * scal / norm_ + 6.0 / rho.powi(2) * a + 1.0 - 0.5 * self.delta)
                * difj
                * difk
                * norm_
            + 16.0 / rho.powi(3) * (dh * norm_.powi(2) + 2.0 * hj * difk)
            + 64.0 / rho.powi(3) * hk * difj
            + 8.0 / rho.powi(2)
                * (16.0 / rho * scal / norm_ - 15.0 / rho.powi(2) * a - 3.0 + 1.5 * self.delta)
                * difj
                * difk
            - 12.0 / rho.powi(2) * (dh * norm_ + hj * difk / norm_)
            - 36.0 / rho.powi(2) * hk * difj / norm_
            + 3.0 / rho
                * (-12.0 / rho * scal / norm_ + 10.0 / rho.powi(2) * a + 3.0 - 1.5 * self.delta)
                * difj
                * difk
                / norm_;
        if the_same {
            dQ_jk = dQ_jk
                + 5.0 * norm_.powi(3) / rho.powi(4)
                    * (-6.0 / rho * scal / norm_ + 6.0 / rho.powi(2) * a + 1.0 - 0.5 * self.delta)
                + 4.0 * norm_.powi(2) / rho.powi(3)
                    * (16.0 / rho * scal / norm_ - 15.0 / rho.powi(2) * a - 3.0 + 1.5 * self.delta)
                + norm_ / rho.powi(2)
                    * (-12.0 / rho * scal / norm_ + 10.0 / rho.powi(2) * a + 3.0
                        - 1.5 * self.delta)
                + self.delta;
        }
        dQ_jk
    }
}

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
