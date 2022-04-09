use std::mem;

pub struct Vec {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec {
    pub fn add(&self, other: &Vec) -> Self {
        Self {
            x: self.x+other.x,
            y: self.y+other.y,
            z: self.z+other.z,
        }
    }

    pub fn sub(&self, other: &Vec) -> Self {
        Self {
            x: self.x-other.x,
            y: self.y-other.y,
            z: self.z-other.z,
        }
    }

    pub fn mul(&self, other: &Vec) -> Self {
        Self {
            x: self.x*other.x,
            y: self.y*other.y,
            z: self.z*other.z,
        }
    }

    pub fn mul_f(&self, d: f64) -> Self {
        Self {
            x: self.x*d,
            y: self.y*d,
            z: self.z*d,
        }
    }

    pub fn norm(&mut self) -> &mut Self {
        let denom = (self.x.powi(2)+self.y.powi(2)+self.z.powi(2)).sqrt().recip();
        self.x *= denom;
        self.y *= denom;
        self.z *= denom;
        self
    }

    pub fn dot(&self, other: &Vec) -> f64 {
        self.x*other.x + self.y*other.y + self.z+other.z
    }

    pub fn modulo(&self, other: &Vec) -> Self {
        Self {
            x: self.y*other.z - self.z*other.y,
            y: self.z*other.x - self.x*other.z,
            z: self.x*other.y - self.y*other.x,
        }
    }
}

pub struct Ray {
    o: Vec,
    d: Vec,
}

impl Ray {
    pub fn new(o: Vec, d: Vec) -> Self {
        Self {
            o,
            d,
        }
    }
}

pub enum Refl {
    Diff,
    Spec,
    Refr,
}

pub struct Sphere {
    rad: f64,
    p: Vec,
    e: Vec,
    c: Vec,
    refl: Refl,
}

impl Sphere {
    pub fn new(rad: f64, p: Vec, e: Vec, c: Vec, refl: Refl) -> Self {
        Self {
            rad,
            p,
            e,
            c,
            refl,
        }
    }

    pub fn intersect(&self, r: &Ray) -> f64 {
        let op = self.p.sub(&r.o);
        let eps = 1e-4f64;
        let b = op.dot(&r.d);
        let mut det = b.powi(2) - op.dot(&op) + self.rad.powi(2);
        if det < 0f64 {
            return 0f64;
        }
        else {
            det = det.sqrt();
        }

        let t = b-det;
        if t > eps {
            return t;
        }
        let t = b+det;
        if t > eps {
            return b+det
        }

        0f64
    }
}

// clamp and toInt functions: use built-ins

fn intersect(spheres: &[Sphere], r: &Ray, t: &mut f64, id: &mut i32) -> bool {
    let inf = 1e20f64;
    *t = inf; // Fix: were originally passed by referece
    for (i, sphere) in spheres.iter().enumerate() {
        let d = sphere.intersect(&r);
        if d != 0f64 && d < *t {
            *t = d;
            *id = i as i32;
        }
    }
    *t < inf
}

fn main() {
}