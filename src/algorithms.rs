#[allow(clippy::too_many_arguments)]
fn rdp_step(x: &[f32], y: &[f32], first: usize, last: usize, epsilon: f32, xnew: &mut [f32], ynew: &mut [f32], count: &mut usize) {
    let mut dmax: f32 = 0.0;
    let mut idx: usize = 0;
    let mut d: f32;
    let m = (y[first] - y[last]) / (x[first] - x[last]);
    let c = y[first] - m * x[first];
    let modulus = (1.0 + m.powi(2)).sqrt();
    for i in (first + 1)..last {
        d = (y[i] - m * x[i] - c).abs() / modulus;
        if d > dmax {
            idx = i;
            dmax = d;
        }
    }
    if dmax > epsilon {
        if idx > first + 1 {
            rdp_step(x, y, first, idx, epsilon, xnew, ynew, count)
        } 
        xnew[*count] = x[idx];
        ynew[*count] = y[idx];
        *count += 1;
        if last > idx + 1 {
            rdp_step(x, y, idx, last, epsilon, xnew, ynew, count)
        } 
    }
}

pub fn ramer_douglas_peucker(x: &[f32], y: &[f32], epsilon: f32) -> Vec<f32> {
    let last = x.len() - 1;
    let mut xnew = [0.0;111];
    let mut ynew = [0.0;111];
    xnew[0] = x[0];
    ynew[0] = y[0];
    let mut count: usize = 1;
    rdp_step(x, y, 0, last, epsilon, &mut xnew, &mut ynew, &mut count);
    xnew[count] = x[last];
    ynew[count] = y[last];
    count += 1;
    [&[count as f32], &xnew[..count], &ynew[..count]].concat()
}