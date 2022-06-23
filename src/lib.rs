mod consts;
mod functions;
mod euler;
mod rk4;
mod other;

#[cfg(test)]
mod tests {
    use crate::euler::euler;
    use crate::rk4::rk4;
    use crate::other::*;

    #[test]
    fn it_works() {
        let soln1 = euler(1.0);
        display_euler(soln1);
        let soln2 = rk4(1.0);
        display_rk4(soln2);
    }
}
