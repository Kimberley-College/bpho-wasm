//* JUST FOR CHECKING VALUES ARE RIGHT */
pub fn _display_rk4(soln: [[f64;111];3]) {
    println!("h       p      t      l");
    for i in 0..23 {
        println!("{:.2}  {:.2}  {:.2}  {:.2}", (i as f64) / 10.0, soln[0][i * 5], soln[1][i * 5], soln[2][i * 5])
    }
}

pub fn _display_euler(soln: [[f64;1101];3]) {
    println!("h       p      t      l");
    for i in 0..23 {
        println!("{:.2}  {:.2}  {:.2}  {:.2}", (i as f64) / 10.0, soln[0][i * 50], soln[1][i * 50], soln[2][i * 50])
    }
}