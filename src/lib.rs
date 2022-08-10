//! A rust implementation of Runge-Kutta's 4th Order and Euler's methods to solve the ODE system laid out in the BPhO Computational Challenge 2022.
extern crate wasm_bindgen;
use wasm_bindgen::prelude::*;

pub mod algorithms;
pub mod schemes;

// initial pressure
const P0: f32 = 1013.25;
// initial temperature
const T0: f32 = 15.0;
// cant do f32 arithmetic in const functions, sooo....
const HVALS: [f32; 111] = [
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.90000004, 1.0, 1.1, 1.2, 1.3000001, 1.4, 1.5,
    1.6, 1.7, 1.8000001, 1.9, 2.0, 2.1000001, 2.2, 2.3, 2.4, 2.5, 2.6000001, 2.7, 2.8, 2.9, 3.0,
    3.1000001, 3.2, 3.3, 3.4, 3.5, 3.6000001, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2000003, 4.3, 4.4, 4.5,
    4.6, 4.7000003, 4.8, 4.9, 5.0, 5.1, 5.2000003, 5.3, 5.4, 5.5, 5.6, 5.7000003, 5.8, 5.9, 6.0,
    6.1, 6.2000003, 6.3, 6.4, 6.5, 6.6, 6.7000003, 6.8, 6.9, 7.0, 7.1, 7.2000003, 7.3, 7.4, 7.5,
    7.6, 7.7000003, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.400001, 8.5, 8.6, 8.7, 8.8, 8.900001, 9.0, 9.1,
    9.2, 9.3, 9.400001, 9.5, 9.6, 9.7, 9.8, 9.900001, 10.0, 10.1, 10.2, 10.3, 10.400001, 10.5,
    10.6, 10.7, 10.8, 10.900001, 11.0,
];

use algorithms::ramer_douglas_peucker;
use schemes::*;

#[wasm_bindgen(js_name = eulerScheme)]
pub fn euler_scheme(u: f32, tolerance: f32) -> Box<[f32]> {
    let soln = euler(u, P0, T0);
    let mut result = Vec::with_capacity(32);
    for i in soln.iter() {
        let v = ramer_douglas_peucker(&HVALS, i, tolerance);
        result.extend(v);
    }
    result.into_boxed_slice()
}

#[wasm_bindgen(js_name = rkScheme)]
pub fn rk4_scheme(u: f32, tolerance: f32) -> Box<[f32]> {
    let soln = rk4(u, P0, T0);
    let mut result = Vec::with_capacity(32);
    for i in soln.iter() {
        let v = ramer_douglas_peucker(&HVALS, i, tolerance);
        result.extend(v);
    }
    result.into_boxed_slice()
}
