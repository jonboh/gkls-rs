use std::ffi::c_double;

use crate::{Domain, GKLSError, Options};

#[allow(dead_code)]
#[link(name = "gkls", kind = "static")]
extern "C" {

    static mut GKLS_dim: usize;
    static mut GKLS_num_minima: usize;
    static mut GKLS_global_dist: f64;
    static mut GKLS_global_radius: f64;
    static mut GKLS_global_value: f64;

    /// allocate boundary vectors
    fn GKLS_domain_alloc() -> i32;
    /// deallocate boundary vectors
    fn GKLS_domain_free();
    /// set default values of the input parameters and allocate the boundary vectors if necessary
    fn GKLS_set_default() -> i32;
    /// deallocate memory needed for the generator
    fn GKLS_free();
    /// test the validity of the input parameters
    fn GKLS_parameters_check() -> i32;
    /// test function generator
    fn GKLS_arg_generate(nf: u32) -> i32;
    /// evaluation of an ND-typed test function
    fn GKLS_ND_func(x: *const f64) -> f64;
    /// evaluation of a D-typed test function
    fn GKLS_D_func(x: *const f64) -> f64;
    /// evaluation of a D2-type test function
    fn GKLS_D2_func(x: *const f64) -> f64;

    /// first order partial derivative of the D-typed test function
    fn GKLS_D_deriv(var_j: usize, x: *const f64) -> f64;
    /// first order partial derivative of the D2-typed test function
    fn GKLS_D2_deriv1(var_j: usize, x: *const f64) -> f64;
    /// second order partial derivative of the D2-typed test function
    fn GKLS_D2_deriv2(var_j: usize, var_k: usize, x: *const f64) -> f64;
    /// gradient of the D-type test function
    fn GKLS_D_gradient(x: *const f64, g: *mut f64) -> i32;
    /// gradient of the D2-type test function fn GKLS_D2_hessian  (double *, double **) -> i32;
    fn GKLS_D2_gradient(x: *const f64, g: *mut f64) -> i32;
    /// Hessian of the D2-type test function
    fn GKLS_D2_hessian(x: *const f64, g: *const *mut f64) -> i32;
}

#[allow(dead_code)]
pub struct CGKLSProblem {
    nf: usize,
    options: Options,
    dim: usize,
    num_minima: usize,
    global_value: f64,
    global_radius: f64,
    global_dist: f64,
}

impl Default for CGKLSProblem {
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

#[allow(dead_code)]
impl CGKLSProblem {
    pub(crate) fn new(
        nf: usize,
        options: Options,
        dim: usize,
        num_minima: usize,
        global_value: f64,
        global_radius: f64,
        global_dist: f64,
    ) -> Result<Self, GKLSError> {
        Ok(CGKLSProblem {
            nf,
            options,
            dim,
            num_minima,
            global_value,
            global_radius,
            global_dist,
        })
    }

    pub(crate) fn nd_func(&self, x: &[f64]) -> f64 {
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap())
        };
        let result = unsafe { GKLS_ND_func(x.as_ptr()) };
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
        result
    }

    pub(crate) fn d_func(&self, x: &[f64]) -> f64 {
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap())
        };
        let result = unsafe { GKLS_D_func(x.as_ptr()) };
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
        result
    }

    /// evaluation of a D2-type test function
    pub(crate) fn d2_func(&self, x: &[f64]) -> f64 {
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap())
        };
        let result = unsafe { GKLS_D2_func(x.as_ptr()) };
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
        result
    }

    /// first order partial derivative of the D-typed test function
    pub(crate) fn d_deriv(&self, var_j: usize, x: &[f64]) -> f64 {
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap())
        };
        let result = unsafe { GKLS_D_deriv(var_j, x.as_ptr()) };
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
        result
    }
    /// first order partial derivative of the D2-typed test function
    pub(crate) fn d2_deriv1(&self, var_j: usize, x: &[f64]) -> f64 {
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap())
        };
        let result = unsafe { GKLS_D2_deriv1(var_j, x.as_ptr()) };
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
        result
    }
    /// second order partial derivative of the D2-typed test function
    pub(crate) fn d2_deriv2(&self, var_j: usize, var_k: usize, x: &[f64]) -> f64 {
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap())
        };
        let result = unsafe { GKLS_D2_deriv2(var_j, var_k, x.as_ptr()) };
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
        result
    }
    /// gradient of the D-type test function
    fn d_gradient(&self, x: &[f64]) -> Vec<f64> {
        let mut g: Vec<f64> = Vec::new();
        g.resize(self.dim, f64::NAN);
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap());
            GKLS_D_gradient(x.as_ptr(), g.as_mut_ptr());
            GKLS_domain_free();
            GKLS_free();
        }
        g
    }
    /// gradient of the D2-type test function fn GKLS_D2_hessian  (double *, double **) -> i32; // Hessian of the D2-type test function
    pub(crate) fn d2_gradient(&self, x: &[f64]) -> Vec<f64> {
        let mut g: Vec<f64> = Vec::new();
        g.resize(self.dim, f64::NAN);
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap());
            GKLS_D2_gradient(x.as_ptr(), g.as_mut_ptr());
            GKLS_domain_free();
            GKLS_free();
        }
        g
    }
    pub(crate) fn d2_hessian(&self, x: &[f64]) -> Vec<Vec<f64>> {
        let mut h: Vec<Vec<f64>> = Vec::new();
        h.resize(self.dim, (0..self.dim).map(|_| f64::NAN).collect());
        let mut row_pointers: Vec<*mut c_double> = Vec::with_capacity(self.dim);
        for row in &mut h {
            row_pointers.push(row.as_mut_ptr());
        }
        unsafe {
            GKLS_dim = self.dim;
            GKLS_num_minima = self.num_minima;
            GKLS_global_dist = self.global_dist;
            GKLS_global_radius = self.global_radius;
            GKLS_global_value = self.global_value;
            GKLS_domain_alloc();
            GKLS_arg_generate(self.nf.try_into().unwrap());
            GKLS_D2_hessian(x.as_ptr(), row_pointers.as_ptr());
            GKLS_domain_free();
            GKLS_free();
        }
        h
    }
}
