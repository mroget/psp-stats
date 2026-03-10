use nalgebra::Matrix3;
use nalgebra::SVD;
use nalgebra::Vector3;
use float_ord::FloatOrd;

fn centroid(points: &[Vector3<f64>]) -> Vector3<f64> {
    points.iter().sum::<Vector3<f64>>() / points.len() as f64
}

fn kabsch_rmsd(p: &[Vector3<f64>], q: &[Vector3<f64>]) -> f64 {
    assert_eq!(p.len(), q.len());

    // 1. Center points
    let c_p = centroid(p);
    let c_q = centroid(q);
    let p_centered: Vec<Vector3<f64>> = p.iter().map(|v| v - c_p).collect();
    let q_centered: Vec<Vector3<f64>> = q.iter().map(|v| v - c_q).collect();

    // 2. Covariance
    let mut cov = Matrix3::<f64>::zeros();
    
    for (pi, qi) in p_centered.iter().zip(q_centered.iter()) {
        cov += pi * qi.transpose();
    }


    // 3. SVD
    let svd = SVD::new(cov, true, true);
    let (u, v_t) = (svd.u.unwrap(), svd.v_t.unwrap());

    // 4. Reflection fix
    let mut d = Matrix3::identity();
    if (u.determinant() * v_t.determinant()) < 0.0 {
        d[(2,2)] = -1.0;
    }

    // 5. Rotation
    let rotation = u * d * v_t;

    // 6. Translation
    let translation = c_p - rotation * c_q;

    // 7. RMSD (after rotation) 
    let mut rmsd_sum = 0.0;
    for (pi, qi) in p.iter().zip(q.iter()) {
        let diff = rotation * qi + translation - pi;
        rmsd_sum += diff.norm_squared();
    }
    let rmsd = (rmsd_sum / p.len() as f64).sqrt();

    rmsd
}


pub fn rmsd(s1 : &[[f64; 3]], s2 : &[[f64; 3]]) -> f64 {
    let p : Vec<Vector3<f64>> = s1.iter().map(|x| Vector3::new(x[0],x[1],x[2])).collect();
    let q : Vec<Vector3<f64>> = s2.iter().map(|x| Vector3::new(x[0],x[1],x[2])).collect();

    kabsch_rmsd(&p,&q)
}

pub fn rmsd_multiple(s1 : &[[f64; 3]], s2 : &Vec<Vec<[f64; 3]>>) -> f64 {
    s2.iter().map(|s| FloatOrd(rmsd(s1,s))).min().unwrap_or(FloatOrd(f64::NAN)).0
}