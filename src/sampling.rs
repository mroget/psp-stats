
use crate::tqdm_;
use pivot_saw::lattice::BaseLattice;
use pivot_saw::lattice::Lattice;
use pivot_saw::walk::SAWIterator;
use pivot_saw::lattice::Tetrahedral;
use crate::dimmerize;
use crate::tqdm;



fn to_f64_coords(sol : &[[i32;3]], scale : f64) -> Vec<[f64;3]> {
    sol.iter().map(|p| p.map(|x| (x as f64)*scale)).collect()
}

fn norm<const D : usize>(p : [i32;D]) -> f64 {
    let mut r = 0.;
    for i in 0..D {
        r = r + (p[i]*p[i]) as f64;
    }
    r.sqrt()
}


pub enum Lat {
    Tetrahedral(f64),
    FCC(f64),
    BCC(f64),
    Cubic(f64),
}


pub enum Method {
    Pivot(usize, usize),
    Dimerize,
    Iterate
}
impl Method {
    fn sample<L: Lattice<i32, 3, N>, const N : usize, F : Fn(Vec<[f64;3]>) -> T, T>(
        self, 
        f : F, 
        len : usize, 
        sample_size : usize, 
        lattice : L, 
        scale : f64, 
        verbose : bool) -> Vec<T> {

        match self {
            Method::Pivot(thermalization_factor, autocorrelation_factor) => {
                let mut pivot =  lattice.get_pivot(len, rand::rng(), thermalization_factor, autocorrelation_factor); 
                tqdm!(0..sample_size, verbose)
                    .map(|_| f(to_f64_coords(&pivot.next().unwrap(), scale)))
                    .collect()
            }
            Method::Dimerize => {
                tqdm!(0..sample_size, verbose)
                    .map(|_| f(to_f64_coords(&dimmerize(&lattice, len), scale)))
                    .collect()
            }
            Method::Iterate => {
                let p0 = [0,0,0];
                let p1 = lattice.neighbors(p0)[0];
                tqdm!(SAWIterator::new(lattice, len, vec![p0,p1]), verbose)
                    .map(|x : Vec<[i32;3]>| f(to_f64_coords(&x, scale)))
                    .collect()
            }
        }
    }
}

pub fn sample_apply<F : Fn(Vec<[f64;3]>) -> T, T>(
        f : F,
        len : usize,
        sample_size : usize,
        method : Method,
        lat : Lat,
        verbose : bool,
    ) -> Vec<T> {
    match lat {
        Lat::Tetrahedral(a) => {
            let lattice = Tetrahedral::new(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(f, len, sample_size, lattice, scale, verbose)
        },
        Lat::Cubic(a) => {
            let lattice = BaseLattice::cubic_grid(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(f, len, sample_size, lattice, scale, verbose)
        },
        Lat::BCC(a) => {
            let lattice = BaseLattice::bcc(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(f, len, sample_size, lattice, scale, verbose)
        },
        Lat::FCC(a) => {
            let lattice = BaseLattice::fcc(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(f, len, sample_size, lattice, scale, verbose)
        },
    }
}
