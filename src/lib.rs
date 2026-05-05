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
mod minimize;

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

fn vec_stats(l : &[f64]) -> (usize,f64,usize) {
    assert!(l.len() != 0);
    let mut argmin = 0;
    let mut argmax = 0;
    let mut avg = l[0];
    for i in 1..l.len() {
        avg += l[i];
        if l[i] < l[argmin] { argmin = i; }
        if l[i] > l[argmax] { argmax = i; }
    }
    avg = avg / (l.len() as f64);
    (argmin, avg, argmax)
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
    ) -> (f64, (f64,f64,f64), (f64,f64,f64), (f64,f64)) {    
    let (e,r) = calculate(seq, gt, cost, sample_size, method, lat, verbose);
    let corr = spearmanr(&e,&r);
    let (amin_e, avg_e, amax_e) = vec_stats(&e);
    let (amin_r, avg_r, amax_r) = vec_stats(&r);
    (corr, (e[amin_e], avg_e, e[amax_e]), (r[amin_r], avg_r, r[amax_r]), (r[amin_e], r[amax_e]))
}


/// A Python module implemented in Rust!
#[pymodule]
mod qpsp_stats {
    use crate::minimize::minimize;
use crate::minimize::Metric;
use crate::Cost;
    use crate::Lat;
    use crate::rmsd_multiple;
    use crate::Method;
    use pyo3::prelude::*;
    use crate::RMSDMultiple;
    
    
    use crate::correlation as rust_correlation;
    use crate::sample_solutions as rust_sample_solutions;
    use crate::stats as rust_stats;


    ///This function calculate several stats for a cost function on a QPSP instance.
    ///
    ///:param seq: The sequence of the protein.
    ///:type seq: str
    ///
    ///:param gt: The groudn truth. A list of 3D strctures (list of lists of 3D vectors).
    ///:type gt: list of floats of shape (k,n,3)
    ///
    ///:param cost: The cost function to be calculated. Can be a string tuple containing lua code and lua function name, a list of weight or an arbitrary python function.
    ///:type cost: (str, str) | dict(dict(str)) | fun(seq : str, coordinates : list(list(f64))) -> f64
    ///
    ///:param sample_size: The number of structures. Ignored if method="iterate".
    ///:type sample_size: int, default: 10000
    ///
    ///:param method: The sampling method.
    ///:type method: str, default: "pivot"
    ///
    ///:param thermalization_factor: The thermalization factor of the pivot algorithm. Ignored for methods other than "pivot".
    ///:type thermalization_factor: int, default: 100
    ///
    ///:param autocorrelation_factor: The autocorrelation factor of the pivot algorithm. Ignored for methods other than "pivot".
    ///:type autocorrelation_factor: int, default: 10
    ///
    ///:param lattice: The 3D lattice used to sample structures.
    ///:type lattice: str, default: "tetrahedral"
    ///
    ///:param arc_length: The length of the lattice's arcs.
    ///:type arc_length: float, default: 3.8
    ///
    ///:param kmin: The minimum sequence distance for interactions. Ignored of cost isn't a dictionnary of energy coefficients.
    ///:type kmin: int, default: 1
    ///
    ///:param dmax: The maximum euclidian distance for interactions. Ignored of cost isn't a dictionnary of energy coefficients.
    ///:type dmax: float, default: 7.8
    ///
    ///:param pbar: If true, show a progress bar.
    ///:type pbar: bool, default: false
    ///
    ///:return: (correlation, (min_cost, avg_cost, max_cost), (min_rmsd, avg_rmsd, max_rmsd), (min_cost_rmsd, max_cost_rmsd))
    ///:rtype: (float, (float, float, float), (float, float, float), (float, float))
    #[pyfunction]
    #[pyo3(signature = (seq, gt, cost, sample_size=10000, method=None, thermalization_factor=100, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, kmin=1, dmax=7.8, pbar=false))]
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
    ) -> PyResult<(f64, (f64,f64,f64), (f64,f64,f64), (f64,f64))> {
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


    ///This function calculate the correlation factor for a cost function on a QPSP instance.
    ///
    ///:param seq: The sequence of the protein.
    ///:type seq: str
    ///
    ///:param gt: The groudn truth. A list of 3D strctures (list of lists of 3D vectors).
    ///:type gt: list of floats of shape (k,n,3)
    ///
    ///:param cost: The cost function to be calculated. Can be a string tuple containing lua code and lua function name, a list of weight or an arbitrary python function.
    ///:type cost: (str, str) | dict(dict(str)) | fun(seq : str, coordinates : list(list(f64))) -> f64
    ///
    ///:param sample_size: The number of structures. Ignored if method="iterate".
    ///:type sample_size: int, default: 10000
    ///
    ///:param method: The sampling method.
    ///:type method: str, default: "pivot"
    ///
    ///:param thermalization_factor: The thermalization factor of the pivot algorithm. Ignored for methods other than "pivot".
    ///:type thermalization_factor: int, default: 100
    ///
    ///:param autocorrelation_factor: The autocorrelation factor of the pivot algorithm. Ignored for methods other than "pivot".
    ///:type autocorrelation_factor: int, default: 10
    ///
    ///:param lattice: The 3D lattice used to sample structures.
    ///:type lattice: str, default: "tetrahedral"
    ///
    ///:param arc_length: The length of the lattice's arcs.
    ///:type arc_length: float, default: 3.8
    ///
    ///:param kmin: The minimum sequence distance for interactions. Ignored of cost isn't a dictionnary of energy coefficients.
    ///:type kmin: int, default: 1
    ///
    ///:param dmax: The maximum euclidian distance for interactions. Ignored of cost isn't a dictionnary of energy coefficients.
    ///:type dmax: float, default: 7.8
    ///
    ///:param pbar: If true, show a progress bar.
    ///:type pbar: bool, default: false
    ///
    ///:return: The correlation coefficient of function "cost" for the protein.
    ///:rtype: float between -1 and 1
    #[pyfunction]
    #[pyo3(signature = (seq, gt, cost, sample_size=10000, method=None, thermalization_factor=100, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, kmin=1, dmax=7.8, pbar=false))]
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
            "enumerate" => {Method::Iterate},
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

    ///This function generates and returns a list of 3D self avoinding walk on a given lattice.
    ///
    ///:param len: The length of the structures.
    ///:type len: int
    ///
    ///:param sample_size: The number of structures. Ignored if method="iterate".
    ///:type sample_size: int, default: 10000
    ///
    ///:param method: The sampling method.
    ///:type method: str, default: "pivot"
    ///
    ///:param thermalization_factor: The thermalization factor of the pivot algorithm. Ignored for methods other than "pivot".
    ///:type thermalization_factor: int, default: 10
    ///
    ///:param autocorrelation_factor: The autocorrelation factor of the pivot algorithm. Ignored for methods other than "pivot".
    ///:type autocorrelation_factor: int, default: 10
    ///
    ///:param lattice: The 3D lattice used to sample structures.
    ///:type lattice: str, default: "tetrahedral"
    ///
    ///:param arc_length: The length of the lattice's arcs.
    ///:type arc_length: float, default: 3.8
    ///
    ///:param pbar: If true, show a progress bar.
    ///:type pbar: bool, default: false
    ///
    ///:return: A list of 3D self avoiding walk on the lattice.
    ///:rtype: list of floats of shape (k,len,3). If the method is not "iterate", then k=sample_size.
    #[pyfunction]
    #[pyo3(signature = (len, sample_size=10000, method="pivot", thermalization_factor=10, autocorrelation_factor=10, lattice="tetrahedral", arc_length=3.8, pbar=false))]
    fn sample_solutions(len : usize,
        sample_size : usize,
        method : &str,
        thermalization_factor : usize, 
        autocorrelation_factor : usize,
        lattice : &str,
        arc_length : f64,
        pbar : bool,
    ) -> PyResult<Vec<Vec<[f64;3]>>> {
            let m = match method {
                "pivot" => {Method::Pivot(thermalization_factor,autocorrelation_factor)},
                "dimerize" => {Method::Dimerize},
                "enumerate" => {Method::Iterate},
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



    ///This function calculates the structure with minimum rmsd aginst ground truth on a given lattice.
    ///
    ///:param gt: The groudn truth. A list of 3D strctures (list of lists of 3D vectors).
    ///:type gt: list of floats of shape (k,n,3)
    ///
    ///:param lattice: The 3D lattice used to sample structures.
    ///:type lattice: str, default: "tetrahedral"
    ///
    ///:param arc_length: The length of the lattice's arcs.
    ///:type arc_length: float, default: 3.8
    ///
    ///:param pbar: If true, show a progress bar.
    ///:type pbar: bool, default: false
    ///
    ///:return: A 3D self avoiding walk on the lattice.
    ///:rtype: list of floats of shape (len,3). **DO NOT RUN THIS FOR LARGE SEQUENCES AS IT ENUMERATES ALL FEASIBLE SOLUTIONS WHICH GROW EXPONENTIALLY WITH THE SEQUENCE LENGTH !!!**
    #[pyfunction]
    #[pyo3(signature = (gt, lattice="tetrahedral", arc_length=3.8, pbar=false))]
    fn minimize_rmsd(gt : Vec<Vec<[f64;3]>>,
        lattice : &str,
        arc_length : f64,
        pbar : bool,
    ) -> PyResult<Vec<[f64;3]>> {
            assert!(gt.len()>0);

            let lat = match lattice {
                    "tetrahedral" => {Lat::Tetrahedral(arc_length)},
                    "cubic" => {Lat::Cubic(arc_length)},
                    "bcc" => {Lat::BCC(arc_length)},
                    "fcc" => {Lat::FCC(arc_length)},
                    x => {panic!("Lattice {:?} is not recognized !", x);}
                };

            Ok(minimize(
                gt[0].len(),
                Metric::RMSD(RMSDMultiple::new(&gt)),
                lat,
                pbar)
            )
        }


    ///This function calculates the structure with minimum cost on a given lattice.
    ///
    ///:param seq: The sequence of the protein.
    ///:type seq: str
    ///
    ///:param cost: The cost function to be calculated. Can be a string tuple containing lua code and lua function name, a list of weight or an arbitrary python function.
    ///:type cost: (str, str) | dict(dict(str)) | fun(seq : str, coordinates : list(list(f64))) -> f64
    ///
    ///:param lattice: The 3D lattice used to sample structures.
    ///:type lattice: str, default: "tetrahedral"
    ///
    ///:param arc_length: The length of the lattice's arcs.
    ///:type arc_length: float, default: 3.8
    ///
    ///:param kmin: The minimum sequence distance for interactions. Ignored of cost isn't a dictionnary of energy coefficients.
    ///:type kmin: int, default: 1
    ///
    ///:param dmax: The maximum euclidian distance for interactions. Ignored of cost isn't a dictionnary of energy coefficients.
    ///:type dmax: float, default: 7.8
    ///
    ///:param pbar: If true, show a progress bar.
    ///:type pbar: bool, default: false
    ///
    ///:return: A 3D self avoiding walk on the lattice.
    ///:rtype: list of floats of shape (len,3). **DO NOT RUN THIS FOR LARGE SEQUENCES AS IT ENUMERATES ALL FEASIBLE SOLUTIONS WHICH GROW EXPONENTIALLY WITH THE SEQUENCE LENGTH !!!**
    #[pyfunction]
    #[pyo3(signature = (seq, cost, lattice="tetrahedral", arc_length=3.8, kmin=1, dmax=7.8, pbar=false))]
    fn minimize_cost(seq : String, 
        cost : &Bound<'_, PyAny>, 
        lattice : &str,
        arc_length : f64,
        kmin : usize,
        dmax : f64,
        pbar : bool,
    ) -> PyResult<Vec<[f64;3]>> {

            let c = Cost::new(cost, seq.clone(), kmin, dmax);
            
            let lat = match lattice {
                    "tetrahedral" => {Lat::Tetrahedral(arc_length)},
                    "cubic" => {Lat::Cubic(arc_length)},
                    "bcc" => {Lat::BCC(arc_length)},
                    "fcc" => {Lat::FCC(arc_length)},
                    x => {panic!("Lattice {:?} is not recognized !", x);}
                };

            Ok(minimize(
                seq.len(),
                Metric::Cost(c),
                lat,
                pbar)
            )
        }

    
    ///A function that calculates the rmsd. This is how the rmsd is calculated when the correlation or other statistics are computed.
    ///
    ///The idea is to get the rmsd between a given structure and a list of structures (usually the models obtained from wetlabs experiments).
    ///
    ///:param sol: A 3D structure (list of 3D vectors).
    ///:type sol: list of floats of shape (n, 3)
    ///
    ///:param gt: A list of 3D strctures (list of lists of 3D vectors).
    ///:type gt: list of floats of shape (k,n,3)
    ///
    ///:return: This function returns the the minimum rmsd between `sol` and each of the structures in `gt`.
    ///:rtype: float
    ///
    ///This function returns the the minimum rmsd between `sol` and each of the structures in `gt`.
    ///This function uses the kabsh algorithm to align the structures.
    ///This function panic if two structures do not have the same length.
    #[pyfunction]
    fn rmsd(sol: Vec<[f64;3]>, gt: Vec<Vec<[f64;3]>>) -> PyResult<f64> {
        let tmp = RMSDMultiple::new(&gt);
        Ok(tmp.calc(&sol))
        //Ok(rmsd_multiple(&sol, &gt))
    }
}
