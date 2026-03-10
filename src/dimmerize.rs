use std::collections::HashSet;
use crate::Lattice;

fn check_collisions_naive<const D : usize>(walk : &[[i32; D]]) -> bool {
    for i in 0..walk.len() {
        for j in (i+1)..walk.len() {
            if walk[i] == walk[j] {
                return true;
            }
        }
    }
    false
}

fn check_collisions<const D : usize>(walk : &[[i32; D]]) -> bool {
    let mut set = HashSet::new();
    for i in 0..walk.len() {
        if set.contains(&walk[i]) {
            return true;
        }
        else {
            set.insert(walk[i]);
        }
    }
    false
}

fn random_saw<const D : usize, const N : usize, T: Lattice<i32, D, N>>(lat : &T, len : usize) -> Vec<[i32; D]> {
    let mut rng = rand::rng();
    let mut ret = lat.random_walk(len, &mut rng);
    while check_collisions_naive(&ret) {
        ret = lat.random_walk(len, &mut rng);
    }
    ret
}

pub fn dimmerize<const D : usize, const N : usize, T: Lattice<i32, D, N>>(lat : &T, len : usize) -> Vec<[i32; D]> {
        if len <= 50 {
            return random_saw(lat, len);
        }

        let n = len;
        let a = n/2;
        let b = n-a+1;

        loop {
            let mut sol1 = dimmerize(lat, a);
            let sol2 = dimmerize(lat, b);
            sol1.extend(&sol2);
            if !check_collisions(&sol1) {
                return sol1;
            }
        }
    }