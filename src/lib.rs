mod consts;
mod functions;
mod euler;
mod rk4;

#[cfg(test)]
mod tests {
    use crate::euler::euler;
    use crate::rk4::rk4;

    #[test]
    fn it_works() {
        let soln1 = euler(1.0);
        println!("{}", soln1[2][1100]);
        let soln2 = rk4(1.0);
        println!("{}", soln2[2][110]);
    }
}
