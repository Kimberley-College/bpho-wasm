//* Important functions */
// calculating relative humidity
pub fn calc_es(t: f32) -> f32 { 6.1121 * ((18.678 - t / 234.5) * (t / (t + 257.14))).exp() }
// calculating dp/dh
pub fn calc_dp(p: f32, ues: f32, tk: f32) -> f32 { -34.171 * (p - 0.37776 * ues) / tk }
// calculating lapse rate
pub fn calc_l(p: f32, ues: f32, tk: f32) -> f32 {
    let r = ues / (p-ues);
    9.7734 * (tk + 5420.3 * r) / (tk * tk + 8400955.5 * r) * tk
}