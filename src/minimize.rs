use crate::Lat;
use crate::sample_apply;
use float_ord::FloatOrd;
use crate::Method;
use crate::RMSDMultiple;
use crate::cost;

pub enum Metric<'a> {
	Cost(cost::Cost<'a>),
	RMSD(RMSDMultiple)
}

impl Metric<'_> {
	fn eval(&self, sol : &Vec<[f64; 3]>) -> f64 {
		match self {
			Metric::Cost(cost) => {
				cost.call(sol)
			},
			Metric::RMSD(rmsd) => {
				rmsd.calc(sol)
			}
		}
	}
}

pub fn minimize(len : usize, metric : Metric<'_>, lat : Lat, verbose : bool) -> Vec<[f64;3]>{
	sample_apply(|x| (x.clone(), FloatOrd(metric.eval(&x))), len, 0, Method::Iterate, lat, verbose)
	.into_iter()
	.min_by_key(|x|x.1)
	.unwrap().0
}