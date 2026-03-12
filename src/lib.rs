use crate::sampling::sample_apply;
use crate::rmsd::RMSDMultiple;
use crate::cost::Cost;
use crate::sampling::Method;
use crate::sampling::Lat;
use crate::dimmerize::dimmerize;
use correlation::spearmanr;
use crate::rmsd::rmsd_multiple;
use pivot_saw::lattice::Lattice;
use pyo3::prelude::*;
use tqdm::tqdm as tqdm_;

mod rmsd;
mod dimmerize;
mod sampling;
mod cost;


macro_rules! tqdm {
    ($iter:expr, $cond:expr) => {
        if $cond {
            itertools::Either::Left(tqdm_($iter))
        } else {
            itertools::Either::Right($iter)
        }
    };
}

pub(crate) use tqdm;

fn vec_stats(l : &[f64]) -> (f64,f64,f64) {
    assert!(l.len() != 0);
    let mut min = l[0];
    let mut max = l[0];
    let mut avg = l[0];
    for i in 1..l.len() {
        avg += l[i];
        if l[i] < min { min = l[i]; }
        if l[i] > max { max = l[i]; }
    }
    avg = avg / (l.len() as f64);
    (min, avg, max)
}


fn calculate(
        seq : String, 
        gt : &Vec<Vec<[f64;3]>>, 
        cost : Cost,
        sample_size : usize,
        method : Method,
        lat : Lat,
        verbose : bool,
    ) -> (Vec<f64>,Vec<f64>) {
    let rmsd_calculator = RMSDMultiple::new(gt);
    let ret = sample_apply(|xyz| 
        {
            (rmsd_calculator.calc(&xyz), cost.call(&xyz))
        }, seq.len(), sample_size, method, lat, verbose);
    let r = ret.iter().map(|x|x.0).collect();
    let e = ret.into_iter().map(|x|x.1).collect();
    (e,r)
}

fn sample_solutions(
        len : usize,
        sample_size : usize,
        method : Method,
        lat : Lat,
        verbose : bool,
    ) -> Vec<Vec<[f64;3]>> {
    sample_apply(|x| x, len, sample_size, method, lat, verbose)
}


fn correlation(
        seq : String, 
        gt : &Vec<Vec<[f64;3]>>, 
        cost : Cost,
        sample_size : usize,
        method : Method,
        lat : Lat,
        verbose : bool,
    ) -> f64{    
    let (e,r) = calculate(seq, gt, cost, sample_size, method, lat, verbose);
    let s = spearmanr(&e,&r);
    s
}


fn stats(
        seq : String, 
        gt : &Vec<Vec<[f64;3]>>, 
        cost : Cost,
        sample_size : usize,
        method : Method,
        lat : Lat,
        verbose : bool,
    ) -> (f64, (f64,f64,f64), (f64,f64,f64)) {    
    let (e,r) = calculate(seq, gt, cost, sample_size, method, lat, verbose);
    let corr = spearmanr(&e,&r);
    let stat_e = vec_stats(&e);
    let stat_r = vec_stats(&r);
    (corr, stat_e, stat_r)
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
    use crate::stats as rust_stats;


    #[pyfunction]
    #[pyo3(signature = (seq, gt, cost, sample_size=10000, method=None, thermalization_factor=10, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, kmin=1, dmax=7.8, pbar=false))]
    fn stats(seq : String, 
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
        pbar : bool,
    ) -> PyResult<(f64, (f64,f64,f64), (f64,f64,f64))> {
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

        Ok(rust_stats(
            seq, 
            &gt, 
            c, 
            sample_size, 
            m,
            lat,
            pbar)
        )
    }



    #[pyfunction]
    #[pyo3(signature = (seq, gt, cost, sample_size=10000, method=None, thermalization_factor=10, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, kmin=1, dmax=7.8, pbar=false))]
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
        pbar : bool,
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
            lat,
            pbar)
        )
    }


    #[pyfunction]
    #[pyo3(signature = (len, sample_size=10000, method=None, thermalization_factor=10, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, pbar=false))]
    fn sample_solutions(len : usize,
        sample_size : usize,
        method : Option<String>,
        thermalization_factor : usize, 
        autocorrelation_factor : usize,
        lattice : &str,
        arc_length : f64,
        pbar : bool,
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
                lat,
                pbar)
            )
        }

    
    /// A function that calculates the rmsd. This is how the rmsd is calculated when the correlation or other statistics are computed.
    /// The idea is to get the rmsd between a given structure and a list of structures (usually the models obtained from wetlabs experiments).
    /// Args:
    ///     - `sol`: A list of 3D vectors (list of length 3 lists).
    ///     - `gt`: The list of structures that acts as ground truth.
    /// This function returns the the minimum rmsd between `sol` and each of the structures in `gt`.
    /// This function uses the kabsh algorithm to allign the structures.
    /// This function panic if two structures do not have the same length.
    ///
    /// Example:
    /// ```python
    /// sol = [[0.,0.,0.], [1.,1.,1.]]
    /// gt = [
    ///        [[0.,0.,0.], [1.,-1.,1.]],
    ///        [[1.,1.,1.], [0.,0.,0.]]
    /// ]
    /// assert(rmsd(sol,gt) <= 1e-5) # rmsd is equal to 0 between sol and gt[1].
    ///```
    #[pyfunction]
    fn rmsd(sol: Vec<[f64;3]>, gt: Vec<Vec<[f64;3]>>) -> PyResult<f64> {
        Ok(rmsd_multiple(&sol, &gt))
    }
}
