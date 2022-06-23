use crate::functions::*;
use crate::consts::*;

pub fn euler(u: f64) -> [[f64; 1101];3] {
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
    [p, t, l]
}