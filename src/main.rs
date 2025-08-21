use rand::prelude::*;
use std::env;

use mlwe::*;

fn main() {
    let mut rng = thread_rng();
    let args: Vec<String> = env::args().collect();

    // encrypt then decrypt each arg as a string
    for arg in &args[1..] {
        println!("input:  {arg}");
        let params = Params {
            n: arg.len() * 8,
            d: 2,
            q: 65535, // Increased q for more noise tolerance
            eta: 1,
        };
        let ((a, t), s) = keygen(&params, &mut rng);

        let (u, v) = encrypt_bytes(&params, &a, &t, arg.as_bytes(), &mut rng);
        let decrypted_bytes = decrypt_bytes(&params, &s, &u, &v);
        let decrypted_arg = String::from_utf8(decrypted_bytes).unwrap();

        println!("output: {}", decrypted_arg);
    }
    
}
