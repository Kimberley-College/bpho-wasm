use super::*;

/// Runge-Kutta 4th Order
#[inline(always)]
pub fn rk4(u: f32, p0: f32, t0: f32) -> [[f32; 111]; 5] {
    // memory allocations
    let mut soln = [[0.0;111];5];
    let (mut t1, mut t2, mut t3, mut t4): (f32, f32, f32, f32);
    let (mut p1, mut p2, mut p3, mut p4): (f32, f32, f32, f32);
    let (mut ues1, mut ues2, mut ues3, mut ues4): (f32, f32, f32, f32);
    let (mut tk1, mut tk2, mut tk3, mut tk4): (f32, f32, f32, f32);
    let (mut kt1, mut kt2, mut kt3, mut kt4): (f32, f32, f32, f32);
    let (mut kp1, mut kp2, mut kp3, mut kp4): (f32, f32, f32, f32);
    // initial conditions
    soln[0][0] = p0;
    soln[1][0] = t0;
    soln[2][0] = calc_l(p0, u * calc_es(t0), t0 + KELVIN);
    soln[3][0] = calc_dew(t0, u);
    soln[4][0] = calc_boil(p0);
    // solving
    for i in 1..111 {
        // first approximaiton
        p1 = soln[0][i - 1];
        t1 = soln[1][i - 1];
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
        soln[0][i] = p1 + (kp1 + 2.0 * kp2 + 2.0 * kp3 + kp4) * DH / 6.0;
        soln[1][i] = t1 - (kt1 + 2.0 * kt2 + 2.0 * kt3 + kt4) * DH / 6.0;
        soln[2][i] = kt1;
        soln[3][i] = calc_dew(soln[1][i], u);
        soln[4][i] = calc_boil(soln[0][i]);
    }
    soln
}

/// Euler's Method
#[inline(always)]
pub fn euler(u: f32, p0: f32, t0: f32) -> [[f32; 111]; 5] {
    let mut soln = [[0.0;111];5];
    let (mut t1, mut ues1, mut tk1, mut tm, mut pm, mut lm): (f32, f32, f32, f32, f32, f32);
    soln[0][0] = p0;
    soln[1][0] = t0;
    soln[2][0] = calc_l(p0, u * calc_es(t0), t0 + KELVIN);
    soln[3][0] = calc_dew(t0, u);
    soln[4][0] = calc_boil(p0);
    for i in 1..111 {
        pm = soln[0][i - 1];
        tm = soln[1][i - 1];
        lm = soln[2][i - 1];
        for _ in 0..10 {
            t1 = tm - lm * DHE;
            ues1 = u * calc_es(t1);
            tk1 = t1 + KELVIN;
            tm = t1;
            pm = pm + DHE * calc_dp(pm, ues1, tk1);
            lm = calc_l(pm, ues1, tk1);
        }
        soln[0][i] = pm;
        soln[1][i] = tm;
        soln[2][i] = lm;
        soln[3][i] = calc_dew(tm, u);
        soln[4][i] = calc_boil(pm);
    }
    soln
}

#[wasm_bindgen(js_name = eulerShort)]
pub fn euler_short(u: f32, p0: f32, t0: f32) -> Box<[f32]> {
    let es = |t: f32| 6.1121 *  ((18.678 - t / 234.5) * (t / (t + 257.14))).exp();
    let dp = |p: f32, ues: f32, tk: f32| -34.171 * (p - 0.37776 * ues) / tk;
    let l = |p: f32, ues: f32, tk: f32| 9.7734 * (tk + 5420.3 * ues / (p - ues)) / (tk * tk + 8400955.5 * ues / (p - ues)) * tk;
    let next = |u: f32, x: [f32; 3]| [x[0] + 0.01 * dp(x[0], u * es(x[1]), x[1] + 273.15), x[1] - 0.01 * x[2], l(x[0], u * es(x[1]), x[1] + 273.15)];
    let mut soln = vec![[p0, t0, l(p0, u * es(t0), t0 + 273.15)]];
    for _ in 1..1101 { soln.push(next(u, *soln.last().unwrap())); }
    soln.concat().into_boxed_slice()
}