use std::fmt::Formatter;
use std::fmt::Display;
use std::collections::HashMap;
use pyo3::prelude::*;
use mlua::prelude::*;

fn dist(a : [f64; 3], b : [f64; 3]) -> f64 {
    (0..2).map(|i| (a[i]-b[i]).powf(2.)).collect::<Vec<f64>>().into_iter().sum::<f64>().sqrt()
}




pub enum Cost<'py> {
    Python(&'py Bound<'py, PyAny>, String),
    Lua(Lua, LuaFunction, String),
    Contact(Vec<f64>, usize, f64)
}
impl Display for Cost<'_> {
    fn fmt(&self, fmt: &mut Formatter<'_>) -> Result<(), std::fmt::Error> { 
        match self {
            Cost::Python(_,_) => {write!(fmt, "a Python function")},
            Cost::Lua(_,_,_) => {write!(fmt, "a Lua function")},
            Cost::Contact(_,_,_) => {write!(fmt, "contact energies")},
        }
        
    }
}

impl<'py> Cost<'py, > {
    pub fn new(cost : &'py pyo3::Bound<'py, pyo3::PyAny>, seq : String, kmin : usize, dmax : f64) -> Cost<'py> {
        match cost.is_callable() {
            false => {
                match cost.extract::<(String,String)>() {
                    Ok((f, code)) => {
                        // Initialize lua env
                        let lua = Lua::new(); 
                        // Run lua code once
                        match lua.load(code).exec() {
                            Ok(_) => {},
                            Err(e) => {panic!("Lua code returned an error when run:\n{}",e)},
                        }
                        // Register the function to call it later
                        let fun: LuaFunction = lua.globals().get(f.clone())
                            .expect(&format!("Lua function {} cannot be found !", f));

                        Cost::Lua(lua, fun, seq)
                    }   
                    Err(_) => {
                        match cost.extract::<HashMap<char,HashMap<char,f64>>>() {
                            Ok(map) => {
                                let l : Vec<char> = seq.chars().collect();
                                let mut ce = vec![];
                                for i in 0..l.len() {
                                    for j in (i+kmin)..l.len() {
                                        ce.push(*map.get(&l[i])
                                            .expect(&format!("Contact energies for amino acid {} could not be found !", l[i]))
                                            .get(&l[j])
                                            .expect(&format!("Contact energy for amino acid's pair {}-{} could not be found !", l[i], l[j])));
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

    pub fn call(&self, xyz : &Vec<[f64;3]>) -> f64 {
        match self {
            Cost::Lua(_lua, cost, seq) => {
                let xyz : &[[f64;3]] = xyz;
                let seq : &str = seq;
                cost.call((seq, xyz))
                .expect("Lua function did not return a floating number !")
            },
            Cost::Python(cost, seq) => {
                cost.call((seq,xyz), None)
                    .expect("Python function failed to run !")
                    .extract::<f64>()
                    .expect("Python function did not return a float !")
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