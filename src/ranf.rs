pub struct Ranf {
    long_lag: usize,
    short_lag: usize,
    separation: usize,
    QUALITY: usize,
    size: usize,
    ran_u: Vec<f64>, // generator state
    counter: usize,
    rnd_num: Vec<f64>, // cached random numbers
}

impl Default for Ranf {
    fn default() -> Self {
        Self::new(100, 37, 70, 1009, 1009, 42)
    }
}

impl Ranf {
    pub fn new(
        long_lag: usize,
        short_lag: usize,
        separation: usize,
        quality: usize,
        size: usize,
        seed: i64,
    ) -> Ranf {
        let ran_u = Self::initialize_state(long_lag, short_lag, separation, size, seed);
        let mut generator = Self {
            long_lag,
            short_lag,
            separation,
            QUALITY: quality,
            size,
            ran_u,
            rnd_num: vec![0.0; size],
            counter: 0,
        };
        generator.ranf_array();
        generator
    }

    fn initialize_state(
        long_lag: usize,
        short_lag: usize,
        separation: usize,
        size: usize,
        seed: i64,
    ) -> Vec<f64> {
        let mut ran_u = vec![0.0; long_lag];
        let mut j: usize;
        let mut t: i64;
        let mut s: i64;
        let mut u = vec![0.0; long_lag + long_lag - 1];
        let mut ul = vec![0.0; long_lag + long_lag - 1];
        let ulp = 2.0f64.powi(-52); // (1.0 / (1 << 30)) / (1 << 22); // 2 to the -52
        let mut ss : f64 = 2.0 * ulp * ((seed & 0x3fffffff) + 2) as f64;
        for j in 0..long_lag {
            u[j] = ss;
            ul[j] = 0.0; // bootstrap the buffer
            ss += ss;
            if ss >= 1.0 {
                ss -= 1.0 - 2.0 * ulp; // cyclic shift of 51 bits
            }
        }
        for j in long_lag..long_lag + long_lag - 1 {
            u[j] = 0.0;
            ul[j] = 0.0;
        }
        u[1] += ulp;
        ul[1] = ulp; // make u[1] (and only u[1]) "odd"
        s = seed & 0x3fffffff;
        t = TryInto::<i64>::try_into(separation).unwrap() - 1;
        while t > 0 {
            for j in (1..long_lag).rev() {
                ul[j + j] = ul[j];
                u[j + j] = u[j]; // "square"
            }
            for j in (long_lag - short_lag + 1..long_lag + long_lag - 1)
                .rev()
                .step_by(2)
            {
                ul[long_lag + long_lag - 1 - j] = 0.0;
                u[long_lag + long_lag - 1 - j] = u[j] - ul[j];
            }
            for j in (long_lag..long_lag + long_lag - 1).rev() {
                if ul[j] != 0.0 {
                    ul[j - (long_lag - short_lag)] = ulp - ul[j - (long_lag - short_lag)];
                    u[j - (long_lag - short_lag)] = mod_sum(u[j - (long_lag - short_lag)], u[j]);
                    ul[j - long_lag] = ulp - ul[j - long_lag];
                    u[j - long_lag] = mod_sum(u[j - long_lag], u[j]);
                }
            }
            if s % 2 == 1 {
                for j in (1..long_lag+1).rev() {
                    ul[j] = ul[j - 1];
                    u[j] = u[j - 1];
                }
                ul[0] = ul[long_lag];
                u[0] = u[long_lag]; // shift the buffer cyclically
                if ul[long_lag] != 0.0 {
                    ul[short_lag] = ulp - ul[short_lag];
                    u[short_lag] = mod_sum(u[short_lag], u[long_lag]);
                }
            }
            if s != 0 {
                s >>= 1;
            } else {
                t -= 1;
            }
        }
        for j in 0..short_lag {
            ran_u[j + long_lag - short_lag] = u[j];
        }
        for j in short_lag..long_lag {
            ran_u[j - short_lag] = u[j];
        }
        ran_u
    }

    fn _ranf_array(&mut self) {
        let mut j = 0;
        for i in 0..self.long_lag {
            self.rnd_num[j] = self.ran_u[i];
            j += 1;
        }
        for i in self.long_lag..self.size {
            self.rnd_num[i] = mod_sum(
                self.rnd_num[i - self.long_lag],
                self.rnd_num[i - self.short_lag],
            );
        }
        let mut i = 0;
        while j < self.size {
            self.ran_u[i] = mod_sum(
                self.rnd_num[j - self.long_lag],
                self.rnd_num[j - self.short_lag],
            );
            i += 1;
            j += 1;
        }
        while i < self.long_lag {
            self.ran_u[i] = mod_sum(
                self.rnd_num[j - self.long_lag],
                self.ran_u[i - self.short_lag],
            );
            i += 1;
            j += 1;
        }
    }

    fn ranf_array(&mut self) {
        let mut j = 0;
        while j < self.long_lag {
            self.rnd_num[j] = self.ran_u[j];
            j += 1;
        }
        while j < self.size {
            self.rnd_num[j] = mod_sum(
                self.rnd_num[j - self.long_lag],
                self.rnd_num[j - self.short_lag],
            );
            j += 1;
        }
        let mut i = 0;
        while i < self.short_lag {
            self.ran_u[i] = mod_sum(
                self.rnd_num[j - self.long_lag],
                self.rnd_num[j - self.short_lag],
            );
            i += 1;
            j += 1;
        }
        while i < self.long_lag {
            self.ran_u[i] = mod_sum(
                self.rnd_num[j - self.long_lag],
                self.ran_u[i - self.short_lag],
            );
            i += 1;
            j += 1;
        }
    }

    pub fn reset(&mut self) {
        self.ranf_array();
        self.counter = 0
    }

    pub fn gen<T>(&mut self) -> f64 {
        if self.counter == self.size {
            self.reset()
        }
        let num = self.rnd_num[self.counter];
        self.counter += 1;
        num
    }
}

fn mod_sum(x: f64, y: f64) -> f64 {
    (x + y) - (x + y).floor()
}

fn is_odd(s: i64) -> bool {
    s & 1 != 0
}
