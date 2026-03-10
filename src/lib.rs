use std::collections::HashMap;
use pyo3::ffi::c_str;
use tqdm::tqdm;
use std::time::Instant;
use pivot_saw::lattice::BaseLattice;
use pivot_saw::walk::SAWIterator;
use crate::dimmerize::dimmerize;
use correlation::spearmanr;
use crate::rmsd::rmsd_multiple;
use pivot_saw::lattice::Lattice;
use pivot_saw::lattice::Tetrahedral;
use pyo3::prelude::*;
use mlua::prelude::*;
use pyo3::types::PyList;

mod rmsd;
mod dimmerize;


fn dist(a : [f64; 3], b : [f64; 3]) -> f64 {
    (0..2).map(|i| (a[i]-b[i]).powf(2.)).collect::<Vec<f64>>().into_iter().sum::<f64>().sqrt()
}


enum Method {
    Pivot(usize, usize),
    Dimerize,
    Iterate
}
impl Method {
    fn sample<L: Lattice<i32, 3, N>, const N : usize>(self, len : usize, sample_size : usize, lattice : L, scale : f64) -> Vec<Vec<[f64;3]>> {
        match self {
            Method::Pivot(thermalization_factor, autocorrelation_factor) => {
                let mut pivot =  lattice.get_pivot(len, rand::rng(), thermalization_factor, autocorrelation_factor); 
                tqdm(0..sample_size).map(|_| to_f64_coords(&pivot.next().unwrap(), scale)).collect()
            }
            Method::Dimerize => {
                let mut ret = vec![vec![[0.;3];len];sample_size];
                for i in 0..sample_size {
                    ret[i] = to_f64_coords(&dimmerize(&lattice, len), scale);
                }
                ret
            }
            Method::Iterate => {
                SAWIterator::new(lattice, len).map(|x| to_f64_coords(&x, scale)).collect()
            }
        }
    }
}

enum Lat {
    Tetrahedral(f64),
    FCC(f64),
    BCC(f64),
    Cubic(f64),
}

enum Cost<'py> {
    Python(&'py Bound<'py, PyAny>, String),
    Lua(Lua, LuaFunction, String),
    Contact(Vec<f64>, usize, f64)
}

impl<'py> Cost<'py, > {
    fn new(cost : &'py pyo3::Bound<'py, pyo3::PyAny>, seq : String, kmin : usize, dmax : f64) -> Cost<'py> {
        match cost.is_callable() {
            false => {
                match cost.extract::<(String,String)>() {
                    Ok((f, code)) => {
                        let lua = Lua::new();
                        lua.load(code).exec().unwrap();
                        let fun: LuaFunction = lua.globals().get(f).unwrap();
                        Cost::Lua(lua, fun, seq)
                    }   
                    Err(_) => {
                        match cost.extract::<HashMap<char,HashMap<char,f64>>>() {
                            Ok(map) => {
                                let l : Vec<char> = seq.chars().collect();
                                let mut ce = vec![];
                                for i in 0..l.len() {
                                    for j in (i+kmin)..l.len() {
                                        ce.push(*map.get(&l[i]).unwrap().get(&l[j]).unwrap())
                                    }
                                }
                                Cost::Contact(ce, kmin, dmax)
                            },
                            Err(_) => {panic!("Cannot recognize {} !", cost)},
                        }
                    },
                }
            },
            true => {
                Cost::Python(cost, seq)
            }
        }
    }

    fn call(&self, xyz : &Vec<[f64;3]>) -> f64 {
        match self {
            Cost::Lua(_lua, cost, seq) => {
                let xyz : &[[f64;3]] = xyz;
                let seq : &str = seq;
                cost.call((seq, xyz)).unwrap()
            },
            Cost::Python(cost, seq) => {
                //println!("{}",  cost.getattr(intern!(cost.py(), "__name__")).unwrap());
                match cost.call((seq,xyz), None).unwrap().extract::<f64>() {
                    Ok(x) => {x},
                    Err(_) => {panic!("Argument cost must return a float !");},
                }
            },
            Cost::Contact(ce, kmin, dmax) => {
                let mut e = 0.;
                let mut k = 0;
                for i in 0..xyz.len() {
                    for j in (i+kmin)..xyz.len() {
                        if dist(xyz[i], xyz[j]) <= *dmax {
                            e+=ce[k];
                        }
                        k+=1;
                    }
                }
                e
            }
        }
    }
}


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


fn sample_solutions(
        len : usize,
        sample_size : usize,
        method : Method,
        lat : Lat,
    ) -> Vec<Vec<[f64;3]>> {
    match lat {
        Lat::Tetrahedral(a) => {
            let lattice = Tetrahedral::new(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(len, sample_size, lattice, scale)
        },
        Lat::Cubic(a) => {
            let lattice = BaseLattice::cubic_grid(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(len, sample_size, lattice, scale)
        },
        Lat::BCC(a) => {
            let lattice = BaseLattice::bcc(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(len, sample_size, lattice, scale)
        },
        Lat::FCC(a) => {
            let lattice = BaseLattice::fcc(1);
            let scale = a/norm(lattice.neighbors([0;3])[0]);
            method.sample(len, sample_size, lattice, scale)
        },
    }
}



fn correlation(
        seq : String, 
        gt : &Vec<Vec<[f64;3]>>, 
        cost : Cost,
        sample_size : usize,
        method : Method,
        lat : Lat,
    ) -> f64{    
    let start = Instant::now();
    let solutions = sample_solutions(seq.len(), sample_size, method, lat);
    println!("sampling: {:?}", start.elapsed());
    
    let start = Instant::now();
    let mut r = vec![];
    let mut e = vec![];
    for xyz in tqdm(solutions.into_iter()) {
        r.push(rmsd_multiple(&xyz, gt));
        e.push(cost.call(&xyz))
    }
    println!("rmsd: {:?}", start.elapsed());
    spearmanr(&e,&r)
}


/// A Python module implemented in Rust.
#[pymodule]
mod qpsp_correlation {
    use crate::Cost;
use crate::Lat;
use crate::rmsd_multiple;
use crate::Method;
use pyo3::prelude::*;
    
    use crate::correlation as rust_correlation;
    use crate::sample_solutions as rust_sample_solutions;

    #[pyfunction]
    #[pyo3(signature = (seq, gt, cost, sample_size=10000, method=None, thermalization_factor=10, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, kmin=1, dmax=7.8))]
    fn correlation(seq : String, 
        gt : Vec<Vec<[f64;3]>>, 
        cost : &Bound<'_, PyAny>, 
        sample_size : usize,
        method : Option<String>,
        thermalization_factor : usize, 
        autocorrelation_factor : usize,
        lattice : &str,
        arc_length : f64,
        kmin : usize,
        dmax : f64,
    ) -> PyResult<f64> {
        let c = Cost::new(cost, seq.clone(), kmin, dmax);
        let m = match method.clone().unwrap_or("pivot".to_string()).as_str() {
            "pivot" => {Method::Pivot(thermalization_factor,autocorrelation_factor)},
            "dimerize" => {Method::Dimerize},
            "iterate" => {Method::Iterate},
            _ => {panic!("Method {:?} is not recognized !", method);}
        };

        let lat = match lattice {
            "tetrahedral" => {Lat::Tetrahedral(arc_length)},
            "cubic" => {Lat::Cubic(arc_length)},
            "bcc" => {Lat::BCC(arc_length)},
            "fcc" => {Lat::FCC(arc_length)},
            x => {panic!("Lattice {:?} is not recognized !", x);}
        };

        Ok(rust_correlation(
            seq, 
            &gt, 
            c, 
            sample_size, 
            m,
            lat)
        )
    }


    #[pyfunction]
    #[pyo3(signature = (len, sample_size=10000, method=None, thermalization_factor=10, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8))]
    fn sample_solutions(len : usize,
        sample_size : usize,
        method : Option<String>,
        thermalization_factor : usize, 
        autocorrelation_factor : usize,
        lattice : &str,
        arc_length : f64,
    ) -> PyResult<Vec<Vec<[f64;3]>>> {
            let m = match method.clone().unwrap_or("pivot".to_string()).as_str() {
                "pivot" => {Method::Pivot(thermalization_factor,autocorrelation_factor)},
                "dimerize" => {Method::Dimerize},
                "iterate" => {Method::Iterate},
                _ => {panic!("Method {:?} is not recognized !", method);}
            };

            let lat = match lattice {
                    "tetrahedral" => {Lat::Tetrahedral(arc_length)},
                    "cubic" => {Lat::Cubic(arc_length)},
                    "bcc" => {Lat::BCC(arc_length)},
                    "fcc" => {Lat::FCC(arc_length)},
                    x => {panic!("Lattice {:?} is not recognized !", x);}
                };

            Ok(rust_sample_solutions(
                len,
                sample_size, 
                m,
                lat)
            )
        }

    /// Formats the sum of two numbers as string.
    #[pyfunction]
    fn rmsd(sol: Vec<[f64;3]>, gt: Vec<Vec<[f64;3]>>) -> PyResult<f64> {
        Ok(rmsd_multiple(&sol, &gt))
    }
}
