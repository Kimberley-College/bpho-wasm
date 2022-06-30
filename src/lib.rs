mod consts;
mod functions;
mod other;
extern crate wasm_bindgen;
use consts::*;
use functions::*;
use wasm_bindgen::prelude::*;

#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

// Runge-Kutta 4th Order
#[wasm_bindgen]
pub fn rk4(u: f64) -> *const [[f64; 111];3] {
    // stack allocations
    let (mut p, mut t, mut l) : ([f64; 111],[f64; 111],[f64; 111]) = ([0.0;111],[0.0;111],[0.0;111]);
    let (mut t1, mut t2, mut t3, mut t4) : (f64, f64, f64, f64);
    let (mut p1, mut p2, mut p3, mut p4) : (f64, f64, f64, f64);
    let (mut ues1, mut ues2, mut ues3, mut ues4) : (f64, f64, f64, f64);
    let (mut tk1, mut tk2, mut tk3, mut tk4) : (f64, f64, f64, f64);
    let (mut kt1, mut kt2, mut kt3, mut kt4) : (f64, f64, f64, f64);
    let (mut kp1, mut kp2, mut kp3, mut kp4) : (f64, f64, f64, f64);
    // initial conditions
    p[0] = P0;
    t[0] = T0;
    l[0] = calc_l(P0, u*calc_es(T0), T0 + KELVIN);
    // solving
    for i in 1..111 {
        // first approximaiton
        p1 = p[i - 1];
        t1 = t[i - 1];
        ues1 = u * calc_es(t1);
        tk1 = t1 + KELVIN;
        kt1 = calc_l(p1, ues1, tk1);
        kp1 = calc_dp(p1, ues1, tk1);
        // second approximation
        t2 = t1 - DH * kt1 / 2.0;
        p2 = p1 + DH * kp1 / 2.0;
        ues2 = u * calc_es(t2);
        tk2 = t2 + KELVIN;
        kt2 = calc_l(p2, ues2, tk2);
        kp2 = calc_dp(p2, ues2, tk2);
        // third approximation
        t3 = t1 - DH * kt2 / 2.0;
        p3 = p1 + DH * kp2 / 2.0;
        ues3 = u * calc_es(t3);
        tk3 = t3 + KELVIN;
        kt3 = calc_l(p3, ues3, tk3);
        kp3 = calc_dp(p3, ues3, tk3);
        // fourth approximation
        t4 = t1 - DH * kt3;
        p4 = p1 + DH * kp3;
        ues4 = u * calc_es(t4);
        tk4 = t4 + KELVIN;
        kt4 = calc_l(p4, ues4, tk4);
        kp4 = calc_dp(p4, ues4, tk4);
        // weighted final approximation
        t[i] = t1 - (kt1 + 2.0 * kt2 + 2.0 * kt3 + kt4) * DH / 6.0;
        p[i] = p1 + (kp1 + 2.0 * kp2 + 2.0 * kp3 + kp4) * DH / 6.0;
        l[i] = kt1;
    }
    let soln: [[f64;111];3] = [p,t,l];
    std::ptr::addr_of!(soln)
 }

 // Euler's Method
 #[wasm_bindgen]
 pub fn euler(u: f64) -> *const [[f64; 1101];3] {
    let (mut p, mut t, mut l) : ([f64; 1101],[f64; 1101],[f64; 1101]) = ([0.0;1101],[0.0;1101],[0.0;1101]);
    let (mut t1, mut ues1, mut tk1) : (f64, f64, f64);
    p[0] = P0;
    t[0] = T0;
    l[0] = calc_l(P0, u * calc_es(T0), T0 + KELVIN);
    for i in 1..1101 {
        t1 = t[i - 1] - l[i - 1] * DHE;
        ues1 = u * calc_es(t1);
        tk1 = t1 + 273.15;
        t[i] = t1;
        p[i] = p[i - 1] + DHE * calc_dp(p[i - 1], ues1, tk1);
        l[i] = calc_l(p[i], ues1, tk1);
    }
    let soln: [[f64;1101];3] = [p,t,l];
    std::ptr::addr_of!(soln)
}
