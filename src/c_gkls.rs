use std::ffi::c_double;

use crate::{Domain, GKLSError, Options};

use self::cbinding::{
    GKLS_D2_deriv1, GKLS_D2_deriv2, GKLS_D2_func, GKLS_D2_gradient, GKLS_D2_hessian, GKLS_D_deriv,
    GKLS_D_func, GKLS_D_gradient, GKLS_ND_func, GKLS_arg_generate, GKLS_domain_alloc,
    GKLS_domain_free, GKLS_free, PointersGKLS, GLOBAL_CVARIABLES,
};

mod cbinding {
    use lazy_static::lazy_static;
    use std::sync::Mutex;

    #[allow(dead_code)]
    #[link(name = "gkls", kind = "static")]
    extern "C" {
        /// allocate boundary vectors
        pub(super) fn GKLS_domain_alloc() -> i32;
        /// deallocate boundary vectors
        pub(super) fn GKLS_domain_free();
        /// set default values of the input parameters and allocate the boundary vectors if necessary
        pub(super) fn GKLS_set_default() -> i32;
        /// deallocate memory needed for the generator
        pub(super) fn GKLS_free();
        /// test the validity of the input parameters
        pub(super) fn GKLS_parameters_check() -> i32;
        /// test function generator
        pub(super) fn GKLS_arg_generate(nf: u32) -> i32;
        /// evaluation of an ND-typed test function
        pub(super) fn GKLS_ND_func(x: *const f64) -> f64;
        /// evaluation of a D-typed test function
        pub(super) fn GKLS_D_func(x: *const f64) -> f64;
        /// evaluation of a D2-type test function
        pub(super) fn GKLS_D2_func(x: *const f64) -> f64;

        /// first order partial derivative of the D-typed test function
        pub(super) fn GKLS_D_deriv(var_j: usize, x: *const f64) -> f64;
        /// first order partial derivative of the D2-typed test function
        pub(super) fn GKLS_D2_deriv1(var_j: usize, x: *const f64) -> f64;
        /// second order partial derivative of the D2-typed test function
        pub(super) fn GKLS_D2_deriv2(var_j: usize, var_k: usize, x: *const f64) -> f64;
        /// gradient of the D-type test function
        pub(super) fn GKLS_D_gradient(x: *const f64, g: *mut f64) -> i32;
        /// gradient of the D2-type test function fn GKLS_D2_hessian  (double *, double **) -> i32;
        pub(super) fn GKLS_D2_gradient(x: *const f64, g: *mut f64) -> i32;
        /// Hessian of the D2-type test function
        pub(super) fn GKLS_D2_hessian(x: *const f64, g: *const *mut f64) -> i32;

        static mut GKLS_dim: usize;
        static mut GKLS_num_minima: usize;
        static mut GKLS_global_dist: f64;
        static mut GKLS_global_radius: f64;
        static mut GKLS_global_value: f64;
    }

    pub(super) struct PointersGKLS {
        pub gkls_dim: &'static mut usize,
        pub gkls_num_minima: &'static mut usize,
        pub gkls_global_dist: &'static mut f64,
        pub gkls_global_radius: &'static mut f64,
        pub gkls_global_value: &'static mut f64,
    }

    // GLOBAL_CVARIABLES Mutex makes sure that the global stare of GKLS is not altered while
    // a computation is being performed. Each computation is expected to allocate and deallocate
    // the global state of GKLS. This is extremely bad performance wise, but we are only interested
    // in this binding for testing purposes. For performance a user should use the Rust
    // implementation.
    lazy_static! {
        pub(super) static ref GLOBAL_CVARIABLES: Mutex<PointersGKLS> = Mutex::new(unsafe {
            PointersGKLS {
                gkls_dim: &mut GKLS_dim,
                gkls_num_minima: &mut GKLS_num_minima,
                gkls_global_dist: &mut GKLS_global_dist,
                gkls_global_radius: &mut GKLS_global_radius,
                gkls_global_value: &mut GKLS_global_value,
            }
        });
    }
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
        .expect("The default problem is valid")
    }
}

#[allow(dead_code)]
impl CGKLSProblem {
    pub(crate) const fn new(
        nf: usize,
        options: Options,
        dim: usize,
        num_minima: usize,
        global_value: f64,
        global_radius: f64,
        global_dist: f64,
    ) -> Result<Self, GKLSError> {
        Ok(Self {
            nf,
            options,
            dim,
            num_minima,
            global_value,
            global_radius,
            global_dist,
        })
    }

    fn allocate_problem(&self, pointers: &mut PointersGKLS) {
        *pointers.gkls_dim = self.dim;
        *pointers.gkls_num_minima = self.num_minima;
        *pointers.gkls_global_dist = self.global_dist;
        *pointers.gkls_global_radius = self.global_radius;
        *pointers.gkls_global_value = self.global_value;
        unsafe {
            GKLS_domain_alloc();
            GKLS_arg_generate(
                self.nf
                    .try_into()
                    .expect("Number functions are expected to be between 1 and 100"),
            )
        };
    }

    // Take pointers as an input to force the need to still have the lock acquired when
    // deallocating the problem
    fn deallocate_problem(_pointers: &PointersGKLS) {
        unsafe {
            GKLS_domain_free();
            GKLS_free();
        }
    }

    pub(crate) fn nd_func(&self, x: &[f64]) -> f64 {
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        let result = unsafe { GKLS_ND_func(x.as_ptr()) };
        Self::deallocate_problem(pointers);
        result
    }

    pub(crate) fn d_func(&self, x: &[f64]) -> f64 {
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        let result = unsafe { GKLS_D_func(x.as_ptr()) };
        Self::deallocate_problem(pointers);
        result
    }

    /// evaluation of a D2-type test function
    pub(crate) fn d2_func(&self, x: &[f64]) -> f64 {
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        let result = unsafe { GKLS_D2_func(x.as_ptr()) };
        Self::deallocate_problem(pointers);
        result
    }

    /// first order partial derivative of the D-typed test function
    pub(crate) fn d_deriv(&self, var_j: usize, x: &[f64]) -> f64 {
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        let result = unsafe { GKLS_D_deriv(var_j, x.as_ptr()) };
        Self::deallocate_problem(pointers);
        result
    }
    /// first order partial derivative of the D2-typed test function
    pub(crate) fn d2_deriv1(&self, var_j: usize, x: &[f64]) -> f64 {
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        let result = unsafe { GKLS_D2_deriv1(var_j, x.as_ptr()) };
        Self::deallocate_problem(pointers);
        result
    }
    /// second order partial derivative of the D2-typed test function
    pub(crate) fn d2_deriv2(&self, var_j: usize, var_k: usize, x: &[f64]) -> f64 {
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        let result = unsafe { GKLS_D2_deriv2(var_j, var_k, x.as_ptr()) };
        Self::deallocate_problem(pointers);
        result
    }
    /// gradient of the D-type test function
    fn d_gradient(&self, x: &[f64]) -> Vec<f64> {
        let mut g: Vec<f64> = Vec::new();
        g.resize(self.dim, f64::NAN);
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        unsafe {
            GKLS_D_gradient(x.as_ptr(), g.as_mut_ptr());
        }
        Self::deallocate_problem(pointers);
        g
    }
    /// gradient of the D2-type test function fn `GKLS_D2_hessian`  (double *, double **) -> i32; // Hessian of the D2-type test function
    pub(crate) fn d2_gradient(&self, x: &[f64]) -> Vec<f64> {
        let mut g: Vec<f64> = Vec::new();
        g.resize(self.dim, f64::NAN);
        let pointers = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        unsafe {
            GKLS_D2_gradient(x.as_ptr(), g.as_mut_ptr());
        }
        Self::deallocate_problem(pointers);
        g
    }
    pub(crate) fn d2_hessian(&self, x: &[f64]) -> Vec<Vec<f64>> {
        let mut h: Vec<Vec<f64>> = Vec::new();
        h.resize(self.dim, (0..self.dim).map(|_| f64::NAN).collect());
        let mut row_pointers: Vec<*mut c_double> = Vec::with_capacity(self.dim);
        for row in &mut h {
            row_pointers.push(row.as_mut_ptr());
        }
        let pointers: &mut PointersGKLS = &mut GLOBAL_CVARIABLES.lock().unwrap();
        self.allocate_problem(pointers);
        unsafe {
            GKLS_D2_hessian(x.as_ptr(), row_pointers.as_ptr());
        }
        Self::deallocate_problem(pointers);
        h
    }
}
