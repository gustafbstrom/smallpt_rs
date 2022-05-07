use std::io::Write;
use std::ops;
use rand;

#[derive(Clone, Copy)]
pub struct Vector {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            x,
            y,
            z,
        }
    }

    // pub fn add(&self, other: &Vector) -> Self {
    //     Self::new(self.x+other.x, self.y+other.y, self.z+other.z)
    // }

    // pub fn sub(&self, other: &Vector) -> Self {
    //     Self::new(self.x-other.x, self.y-other.y, self.z-other.z)
    // }

    // pub fn mul(&self, d: f64) -> Self {
    //     Self::new(self.x*d, self.y*d, self.z*d)
    // }

    // pub fn mult(&self, other: &Vector) -> Self {
    //     Self::new(self.x*other.x, self.y*other.y, self.z*other.z)
    // }

    pub fn norm(&mut self) -> Vector {
        let denom = (self.x.powi(2)+self.y.powi(2)+self.z.powi(2)).sqrt().recip();
        self.x *= denom;
        self.y *= denom;
        self.z *= denom;
        *self
    }

    pub fn dot(&self, other: &Vector) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }

    // pub fn modulo(&self, other: &Vector) -> Self {
    //     Self {
    //         x: self.y*other.z - self.z*other.y,
    //         y: self.z*other.x - self.x*other.z,
    //         z: self.x*other.y - self.y*other.x,
    //     }
    // }
}

impl ops::Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, rhs: Vector) -> Self::Output {
        Self::new(self.x+rhs.x, self.y+rhs.y, self.z+rhs.z)
    }
}

impl ops::Sub<Vector> for Vector {
    type Output = Vector;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self::new(self.x-rhs.x, self.y-rhs.y, self.z-rhs.z)
    }
}

impl ops::Mul<Vector> for Vector {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        Self::new(self.x*rhs.x, self.y*rhs.y, self.z*rhs.z)
    }
}

impl ops::Mul<f64> for Vector {
    type Output = Vector;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.x*rhs, self.y*rhs, self.z*rhs)
    }
}

impl ops::Rem<Vector> for Vector {
    type Output = Vector;

    fn rem(self, rhs: Vector) -> Self::Output {
        Self {
            x: self.y*rhs.z - self.z*rhs.y,
            y: self.z*rhs.x - self.x*rhs.z,
            z: self.x*rhs.y - self.y*rhs.x,
        }
    }
}

impl ops::Neg for Vector {
    type Output = Vector;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

#[derive(Clone, Copy)]
pub struct Ray {
    o: Vector,
    d: Vector,
}

impl Ray {
    pub fn new(o: Vector, d: Vector) -> Self {
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
    p: Vector,
    e: Vector,
    c: Vector,
    refl: Refl,
}

impl Sphere {
    pub fn new(rad: f64, p: Vector, e: Vector, c: Vector, refl: Refl) -> Self {
        Self {
            rad,
            p,
            e,
            c,
            refl,
        }
    }

    pub fn intersect(&self, r: &Ray) -> f64 {
        let op = self.p - r.o;
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

fn intersect(spheres: &[Sphere], r: &Ray, t: &mut f64, id: &mut usize) -> bool {
    const INF: f64 = 1e20f64;
    *t = INF;
    for (i, sphere) in spheres.iter().enumerate() {
        let d = sphere.intersect(&r);
        if d != 0f64 && d < *t {
            *t = d;
            *id = i;
        }
    }
    *t < INF
}

fn radiance(spheres: &[Sphere], r: &Ray, depth: i32) -> Vector {
    let mut t = 0f64;
    let mut id = 0usize;
    if !intersect(&spheres, r, &mut t, &mut id) {
        return Vector::new(0.0, 0.0, 0.0);
    }
    let ref obj = spheres[id];
    let x = r.o + r.d*t;
    let n = (x-obj.p).norm();
    let nl = if n.dot(&r.d) < 0.0 { n } else { -n };
    let mut f = obj.c;
    let p = f.x.max(f.y).max(f.z); // max refl
    
    let depth = depth + 1;
    if depth > 5 {
        if rand::random::<f64>() < p {
            f = f*(p.recip());
        }
        else {
            return obj.e;
        }
    }

    match obj.refl {
        Refl::Diff => {
            let r1 = 2.0*std::f64::consts::PI*rand::random::<f64>();
            let r2 = rand::random::<f64>();
            let r2s = r2.sqrt();
            let w = nl;
            let u = if w.x.abs() > 0.1 {
                Vector::new(0., 1., 0.)
            }
            else {
                Vector::new(1., 0., 0.)
            };
            let mut u = (u % w).norm();
            let v = w % u;
            let d = ((u*r1.cos() + v*r1.sin())*r2s + w*(1.-r2).sqrt()).norm();
            return obj.e + (f*radiance(spheres, &Ray::new(x, d), depth));
        },
        Refl::Spec => {
            return obj.e + f*radiance(&spheres, &Ray::new(x, r.d - n*2.0*n.dot(&r.d)), depth);
        }
        _ => ()
    }
    
    let refl_ray = Ray::new(x, r.d - n*2.0*n.dot(&r.d));
    let into = n.dot(&nl) > 0f64;
    let nc = 1f64;
    let nt = 1.5f64;
    let nnt = if into {nc/nt} else {nt/nc};
    let ddn = r.d.dot(&nl);
    let cos2t = 1.0 - nnt.powi(2) * (1.0-ddn.powi(2));
    if cos2t < 0f64 {
        return obj.e + f*radiance(spheres,&refl_ray, depth);
    }

    let direction = if into {1f64} else {-1f64};
    let tdir = (r.d*nnt - n*(direction*(ddn*nnt + cos2t.sqrt()))).norm();
    let a = nt-nc;
    let b = nt+nc;
    let r0 = a.powi(2)/b.powi(2);
    let c = 1f64 - (if into {-ddn} else {tdir.dot(&n)});
    let re = r0 + (1f64-r0)*c.powi(5);
    let tr = 1f64 - re;
    let p = 0.25f64 + 0.5*re;
    let rp = re/p;
    let tp = tr/(1f64-p);

    let res = if depth > 2 {
        if rand::random::<f64>() < p {
            radiance(&spheres, &refl_ray, depth)*rp
        }
        else {
            radiance(&spheres, &Ray::new(x, tdir), depth)*tp
        }
    }
    else {
        radiance(&spheres, &refl_ray, depth)*re + radiance(&spheres, &Ray::new(x, tdir), depth)*tr
    };

    obj.e + f*res
}

fn write_to_file(c: &Vec<Vector>, w: usize, h: usize, path: &str) {
    let mut f = std::fs::File::create(path).expect("Could not open file");
    f.write_all(format!("P3\n{} {}\n{}\n", w, h, 255).as_bytes()).unwrap();
    c.iter().for_each(|v| {
        let func = |x: f64| {(x.clamp(0.0, 1.0).powf(1.0/2.2)*255.0 + 0.5) as u8};
        let x = func(v.x);
        let y = func(v.y);
        let z = func(v.z);
        f.write_all(format!("{} {} {} ", x, y ,z).as_bytes()).unwrap();
    });
}

pub fn main() {
    // Scene: radius, position, emission, color, material
    let spheres = [
        Sphere::new(
            1e5, 
            Vector::new(1e5+1.0,40.8,81.6), 
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(0.75, 0.25, 0.25), 
            Refl::Diff
        ), //Left
        Sphere::new(
            1e5,
            Vector::new(-1e5+99.0,40.8,81.6),
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.25,0.25,0.75), 
            Refl::Diff
        ), //Right
        Sphere::new(
            1e5,
            Vector::new(50.0,40.8, 1e5), 
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.75,0.75,0.75), 
            Refl::Diff
        ), //Back
        Sphere::new(
            1e5,
            Vector::new(50.0, 40.8, -1e5+170.0), 
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.0, 0.0, 0.0), 
            Refl::Diff
        ), //Front
        Sphere::new(
            1e5,
            Vector::new(50.0, 1e5, 81.6), 
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.75, 0.75, 0.75), 
            Refl::Diff
        ), //Bottom
        Sphere::new(
            1e5,
            Vector::new(50.0, -1e5+81.6, 81.6), 
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.75,0.75,0.75), 
            Refl::Diff
        ), //Top
        Sphere::new(
            16.5,
            Vector::new(27.0,16.5,47.0), 
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.999, 0.999, 0.999), 
            Refl::Spec
        ), //Mirror
        Sphere::new(
            600.0,
            Vector::new(50.0, 681.6-0.27, 81.6), 
            Vector::new(12.0, 12.0, 12.0), 
            Vector::new(0.0, 0.0, 0.0), 
            Refl::Diff
        ), //Lite
        Sphere::new(
            16.5,
            Vector::new(73.0 ,16.5, 78.0), 
            Vector::new(0.0, 0.0, 0.0), 
            Vector::new(0.999,0.999,0.999), 
            Refl::Refr
        ), //Glass
    ];
    
    const W: usize = 1024;
    const H: usize = 768;
    let args: Vec<String> = std::env::args().collect();
    let samps = if args.len() == 2 {
        args[1].parse::<i32>().unwrap()/4
    }
    else {
        1
    };
    let cam = Ray::new(Vector::new(50., 52., 295.6), Vector::new(0f64, -0.042612, -1f64).norm());
    let cx = Vector::new((W as f64)*0.5135/(H as f64), 0., 0.);
    let cy = (cx%cam.d).norm()*0.5135;
    let mut c = vec![Vector::new(0.0, 0.0, 0.0); W*H];
    
    // for loop here
    for y in 0..H {
        eprint!("\rRendering ({} spp) {:5.2}%", samps*4, 100.*(y as f32)/((H as f32)-1.));
        for x in 0..W {
            let i = (H-y-1)*W+x;
            for sy in 0..2 {
                for sx in 0..2 {
                    let mut r = Vector::new(0., 0., 0.);
                    for _ in 0..samps {
                        let r1 = 2. * rand::random::<f64>();
                        let dx = if r1 < 1.0 {
                            r1.sqrt()-1.0
                        }
                        else {
                            1.0 - (2.0-r1).sqrt()
                        };
                        let r2 = 2. * rand::random::<f64>();
                        let dy = if r2 < 1.0 {
                            r2.sqrt()-1.0
                        }
                        else {
                            1.0 - (2.0-r2).sqrt()
                        };
                        let mut d =
                            cx*(((sx as f64+0.5+dx)/2.0 + x as f64)/W as f64 - 0.5)
                            + cy*(((sy as f64+0.5+dy)/2.0 + y as f64)/H as f64 - 0.5)
                            + cam.d;
                        r = r
                            + radiance(&spheres,&Ray::new(cam.o + d*140., d.norm()), 0) * (samps as f64).recip();
                    }
                    let v = Vector::new(
                        r.x.clamp(0.0, 1.0), 
                        r.y.clamp(0.0, 1.0),
                        r.z.clamp(0.0, 1.0),
                    )*0.25;
                    c[i] = c[i]+ v;
                }
            }
        }
    }
    let path = "image.ppm";
    println!("\nWriting to file {}", path);
    write_to_file(&c, W, H, path);
}
