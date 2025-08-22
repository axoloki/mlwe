use bitvec::prelude::*;
use rand::prelude::*;
//use std::ops::{Add, Neg, Mul, Rem};
use std::ops::{Add, BitOr, Mul, Neg, Sub};

/// Modular reduction to [0, q).
fn modq(x: i64, q: i64) -> i64 {
    ((x % q) + q) % q
}

/// Minimum distance between a and b modulo q.
fn min_dist(a: i64, b: i64, q: i64) -> i64 {
    let diff = modq(a - b, q);
    diff.min(q - diff)
}

#[derive(Clone, Debug, PartialEq)]
/// A polynomial with integer coefficients
pub struct Polynomial {
    pub coeffs: Vec<i64>,
    n: usize,
    q: i64,
}

impl Polynomial {
    fn new(n: usize, q: i64) -> Self {
        Polynomial {
            coeffs: vec![0i64; n],
            n,
            q,
        }
    }

    fn add_impl(&self, other: &Self) -> Self {
        let mut res = Self::new(self.n, self.q);
        for i in 0..self.n {
            res.coeffs[i] = modq(self.coeffs[i] + other.coeffs[i], self.q);
        }
        res
    }

    fn sub_impl(&self, other: &Self) -> Self {
        let mut res = Self::new(self.n, self.q);
        for i in 0..self.n {
            res.coeffs[i] = modq(self.coeffs[i] - other.coeffs[i], self.q);
        }
        res
    }

    fn neg_impl(&self) -> Self {
        let mut res = Self::new(self.n, self.q);
        for i in 0..self.n {
            res.coeffs[i] = modq(-self.coeffs[i], self.q);
        }
        res
    }

    fn mul_impl(&self, other: &Self) -> Self {
        let mut res_coeffs = vec![0i64; 2 * self.n - 1];
        for i in 0..self.n {
            for j in 0..self.n {
                res_coeffs[i + j] =
                    modq(res_coeffs[i + j] + self.coeffs[i] * other.coeffs[j], self.q);
            }
        }
        // Reduce mod x^n + 1
        for i in (self.n..(2 * self.n - 1)).rev() {
            res_coeffs[i - self.n] = modq(res_coeffs[i - self.n] - res_coeffs[i], self.q);
        }
        let mut res = Self::new(self.n, self.q);
        for i in 0..self.n {
            res.coeffs[i] = modq(res_coeffs[i], self.q);
        }
        res
    }
}

impl Add<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn add(self, other: &Polynomial) -> Self::Output {
        Polynomial::add_impl(self, other)
    }
}

impl Add<Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn add(self, other: Polynomial) -> Self::Output {
        Polynomial::add_impl(self, &other)
    }
}

impl Add<&Polynomial> for Polynomial {
    type Output = Polynomial;

    fn add(self, other: &Polynomial) -> Self::Output {
        Polynomial::add_impl(&self, other)
    }
}

impl Add<Polynomial> for Polynomial {
    type Output = Polynomial;

    fn add(self, other: Polynomial) -> Self::Output {
        Polynomial::add_impl(&self, &other)
    }
}

impl Sub<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn sub(self, other: &Polynomial) -> Self::Output {
        Polynomial::sub_impl(self, other)
    }
}

impl Sub<Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn sub(self, other: Polynomial) -> Self::Output {
        Polynomial::sub_impl(self, &other)
    }
}

impl Sub<&Polynomial> for Polynomial {
    type Output = Polynomial;

    fn sub(self, other: &Polynomial) -> Self::Output {
        Polynomial::sub_impl(&self, other)
    }
}

impl Sub<Polynomial> for Polynomial {
    type Output = Polynomial;

    fn sub(self, other: Polynomial) -> Self::Output {
        Polynomial::sub_impl(&self, &other)
    }
}

impl Neg for &Polynomial {
    type Output = Polynomial;

    fn neg(self) -> Self::Output {
        Polynomial::neg_impl(self)
    }
}

impl Neg for Polynomial {
    type Output = Polynomial;

    fn neg(self) -> Self::Output {
        Polynomial::neg_impl(&self)
    }
}

impl Mul<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn mul(self, other: &Polynomial) -> Self::Output {
        Polynomial::mul_impl(self, other)
    }
}

#[derive(Clone)]
/// A vector of polynomials
pub struct Vector {
    polys: Vec<Polynomial>,
    d: usize,
    n: usize,
    q: i64,
}

impl Vector {
    fn new(d: usize, n: usize, q: i64) -> Self {
        Vector {
            polys: vec![Polynomial::new(n, q); d],
            d,
            n,
            q,
        }
    }

    fn add_impl(&self, other: &Self) -> Self {
        let mut res = Self::new(self.d, self.n, self.q);
        for i in 0..self.d {
            res.polys[i] = &self.polys[i] + &other.polys[i];
        }
        res
    }

    fn sub_impl(&self, other: &Self) -> Self {
        let mut res = Self::new(self.d, self.n, self.q);
        for i in 0..self.d {
            res.polys[i] = &self.polys[i] - &other.polys[i];
        }
        res
    }

    fn inner(&self, other: &Self) -> Polynomial {
        let mut res = Polynomial::new(self.n, self.q);
        for i in 0..self.d {
            res = res + (&self.polys[i] * &other.polys[i]);
        }
        res
    }
}

impl Add<&Vector> for &Vector {
    type Output = Vector;

    fn add(self, other: &Vector) -> Self::Output {
        Vector::add_impl(self, other)
    }
}

impl Add<Vector> for &Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Self::Output {
        Vector::add_impl(self, &other)
    }
}

impl Add<&Vector> for Vector {
    type Output = Vector;

    fn add(self, other: &Vector) -> Self::Output {
        Vector::add_impl(&self, other)
    }
}

impl Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Self::Output {
        Vector::add_impl(&self, &other)
    }
}

impl Sub<&Vector> for &Vector {
    type Output = Vector;

    fn sub(self, other: &Vector) -> Self::Output {
        Vector::sub_impl(self, other)
    }
}

impl Sub<Vector> for &Vector {
    type Output = Vector;

    fn sub(self, other: Vector) -> Self::Output {
        Vector::sub_impl(self, &other)
    }
}

impl Sub<&Vector> for Vector {
    type Output = Vector;

    fn sub(self, other: &Vector) -> Self::Output {
        Vector::sub_impl(&self, other)
    }
}

impl Sub<Vector> for Vector {
    type Output = Vector;

    fn sub(self, other: Vector) -> Self::Output {
        Vector::sub_impl(&self, &other)
    }
}

impl BitOr<&Vector> for &Vector {
    type Output = Polynomial;

    fn bitor(self, other: &Vector) -> Self::Output {
        Vector::inner(self, other)
    }
}

impl BitOr<Vector> for &Vector {
    type Output = Polynomial;

    fn bitor(self, other: Vector) -> Self::Output {
        Vector::inner(self, &other)
    }
}

impl BitOr<&Vector> for Vector {
    type Output = Polynomial;

    fn bitor(self, other: &Vector) -> Self::Output {
        Vector::inner(&self, other)
    }
}

impl BitOr<Vector> for Vector {
    type Output = Polynomial;

    fn bitor(self, other: Vector) -> Self::Output {
        Vector::inner(&self, &other)
    }
}

#[derive(Clone)]
/// A matrix of polynomials
pub struct Matrix {
    rows: Vec<Vec<Polynomial>>,
    d: usize,
    n: usize,
    q: i64,
}

impl Matrix {
    fn new(d: usize, n: usize, q: i64) -> Self {
        Matrix {
            rows: vec![vec![Polynomial::new(n, q); d]; d],
            d,
            n,
            q,
        }
    }

    /// Computes A * v.
    fn matmul(&self, v: &Vector) -> Vector {
        let mut res = Vector::new(self.d, self.n, self.q);
        for i in 0..self.d {
            let mut sum = Polynomial::new(self.n, self.q);
            for j in 0..self.d {
                sum = sum + (&self.rows[i][j] * &v.polys[j]);
            }
            res.polys[i] = sum;
        }
        res
    }

    /// Computes A^T * v.
    fn matmul_transposed(&self, v: &Vector) -> Vector {
        let mut res = Vector::new(self.d, self.n, self.q);
        for j in 0..self.d {
            let mut sum = Polynomial::new(self.n, self.q);
            for i in 0..self.d {
                sum = sum + (&self.rows[i][j] * &v.polys[i]);
            }
            res.polys[j] = sum;
        }
        res
    }
}

impl Mul<&Vector> for &Matrix {
    type Output = Vector;

    fn mul(self, other: &Vector) -> Self::Output {
        Matrix::matmul(self, other)
    }
}

impl Mul<Vector> for &Matrix {
    type Output = Vector;

    fn mul(self, other: Vector) -> Self::Output {
        Matrix::matmul(self, &other)
    }
}

impl Mul<&Vector> for Matrix {
    type Output = Vector;

    fn mul(self, other: &Vector) -> Self::Output {
        Matrix::matmul(&self, other)
    }
}

impl Mul<Vector> for Matrix {
    type Output = Vector;

    fn mul(self, other: Vector) -> Self::Output {
        Matrix::matmul(&self, &other)
    }
}

/// Parameters for M-LWE cryptography
pub struct Params {
    /// Polynomial degree (must be power of 2 for x^n + 1).
    pub n: usize,
    /// Module rank.
    pub d: usize,
    /// Modulus.
    pub q: i64,
    /// Eta for ternary sampling (-eta to eta).
    pub eta: i64,
}

fn sample_uniform_poly<R: Rng>(rng: &mut R, n: usize, q: i64) -> Polynomial {
    let mut p = Polynomial::new(n, q);
    for i in 0..n {
        p.coeffs[i] = rng.gen_range(0..q);
    }
    p
}

fn sample_small_poly<R: Rng>(rng: &mut R, n: usize, eta: i64, q: i64) -> Polynomial {
    let mut p = Polynomial::new(n, q);
    for i in 0..n {
        let c: i64 = rng.gen_range(-eta..=eta);
        p.coeffs[i] = modq(c, q);
    }
    p
}

pub fn sample_binary_poly<R: Rng>(rng: &mut R, n: usize, q: i64) -> Polynomial {
    let mut p = Polynomial::new(n, q);
    for i in 0..n {
        p.coeffs[i] = rng.gen_range(0..=1);
    }
    p
}

fn sample_uniform_matrix<R: Rng>(rng: &mut R, d: usize, n: usize, q: i64) -> Matrix {
    let mut mat = Matrix::new(d, n, q);
    for i in 0..d {
        for j in 0..d {
            mat.rows[i][j] = sample_uniform_poly(rng, n, q);
        }
    }
    mat
}

fn sample_small_vec<R: Rng>(rng: &mut R, d: usize, n: usize, eta: i64, q: i64) -> Vector {
    let mut v = Vector::new(d, n, q);
    for i in 0..d {
        v.polys[i] = sample_small_poly(rng, n, eta, q);
    }
    v
}

/// Generates public key (A, t) and secret key s.
pub fn keygen<R: Rng>(params: &Params, rng: &mut R) -> ((Matrix, Vector), Vector) {
    let a = sample_uniform_matrix(rng, params.d, params.n, params.q);
    let s = sample_small_vec(rng, params.d, params.n, params.eta, params.q);
    let e = sample_small_vec(rng, params.d, params.n, params.eta, params.q);
    let t = (&a * &s) + &e;
    ((a, t), s)
}

/// Encrypts a polynomial message m (with small coefficients) using public key (A, t).
pub fn encrypt<R: Rng>(
    params: &Params,
    a: &Matrix,
    t: &Vector,
    m: &Polynomial,
    rng: &mut R,
) -> (Vector, Polynomial) {
    let r = sample_small_vec(rng, params.d, params.n, params.eta, params.q);
    let e1 = sample_small_vec(rng, params.d, params.n, params.eta, params.q);
    // u = A^T*r + e1 (mod q)
    let u = a.matmul_transposed(&r) + &e1;
    let e2 = sample_small_poly(rng, params.n, params.eta, params.q);
    // Scale message by floor(q/2) to encode
    let mut m_scaled = Polynomial::new(params.n, params.q);
    let scale = params.q / 2;
    for i in 0..params.n {
        m_scaled.coeffs[i] = modq(m.coeffs[i] * scale, params.q);
    }
    // v = <t, r> + m (mod q) [scaled to q/2]
    let v = ((t | r) + &e2) + m_scaled;
    (u, v)
}

/// Decrypts the ciphertext (u, v) using secret key s. Returns the message polynomial.
pub fn decrypt(params: &Params, s: &Vector, u: &Vector, v: &Polynomial) -> Polynomial {
    // m ~= v - <u, s> (mod q)
    // This works, because
    // v - <u, s> = <t, r>           + e2 + m - <u, s>
    //            = <(A * s) + e, r> + e2 + m - <A^T * r + e1, s>
    //            = <(A * s) + e, r> + e2 + m - <A^T * r + e1, s>
    //           ~= <(A * s), r>     + e2 + m - <A^T * r, s>
    //           ~= ((A * s)^T * r)  + e2 + m - ((A^T * r)^T * s)
    //           ~= ((A * s)^T * r)  + e2 + m - ((A * r) * s)
    //           ~= ((A * s)^T * r)  + e2 + m - ((A * r^T) * s)
    //           ~= ((A * s)^T * r)  + e2 + m - (r^T * (A * s))
    //           ~= ((A * s)^T * r)  + e2 + m - ((A * s) * r)
    //           ~=                    e2 + m
    let approx = v - (u | s);
    // Decode each coefficient by rounding to nearest multiple of floor(q/2)
    let mut m = Polynomial::new(params.n, params.q);
    let q2 = params.q / 2;
    for i in 0..params.n {
        let c = approx.coeffs[i];
        let dist_to_0 = min_dist(c, 0, params.q);
        let dist_to_q2 = min_dist(c, q2, params.q);
        m.coeffs[i] = if dist_to_q2 < dist_to_0 { 1 } else { 0 };
    }
    m
}

/// Encrypts bytes using public key (A, t).  Assumes (params.n == bytes.len() * 8)
pub fn encrypt_bytes<R: Rng>(
    params: &Params,
    a: &Matrix,
    t: &Vector,
    bytes: &[u8],
    rng: &mut R,
) -> (Vector, Polynomial) {
    //println!("bytes.len {}", bytes.len());
    let mut message = sample_binary_poly(rng, params.n, params.q);
    //println!("message.len {}", message.coeffs.len());

    for i in 0..bytes.len() {
        let bits = bytes[i].view_bits::<Lsb0>();
        //println!("byte {} bits {}", bytes[i], bits);
        for (j, bit) in bits.iter().enumerate() {
            //println!("bit {j}: {bit}");
            message.coeffs[8 * i + j] = if *bit { 1 } else { 0 };
        }
    }

    encrypt(params, a, t, &message, rng)
}

/// Decrypts the ciphertext (u, v) using secret key s. Returns the polynomial converted to bytes.
pub fn decrypt_bytes(params: &Params, s: &Vector, u: &Vector, v: &Polynomial) -> Vec<u8> {
    let decrypted = decrypt(params, s, u, v);
    //println!("decrypted.len {}", decrypted.coeffs.len());

    let mut bytes = Vec::new();
    bytes.resize(decrypted.coeffs.len() / 8, 0u8);

    for i in 0..bytes.len() {
        let bits = bytes[i].view_bits_mut::<Lsb0>();
        for (j, mut bit) in bits.iter_mut().enumerate() {
            *bit = decrypted.coeffs[8 * i + j] != 0;
        }
    }

    bytes
}

/// Test function to verify encryption/decryption correctness over multiple trials.
pub fn test_encryption<R: Rng>(params: &Params, trials: usize, rng: &mut R) -> usize {
    let mut failures = 0;
    let ((a, t), s) = keygen(params, rng);
    for _ in 0..trials {
        let m = sample_binary_poly(rng, params.n, params.q);
        let (u, v) = encrypt(params, &a, &t, &m, rng);
        let decrypted = decrypt(params, &s, &u, &v);
        if m.coeffs != decrypted.coeffs {
            failures += 1;
            println!(
                "Failure: message = {:?}, decrypted = {:?}",
                m.coeffs, decrypted.coeffs
            );
        }
    }
    failures
}

/// Test function to verify encryption/decryption correctness over multiple trials.
pub fn test_encrypt_bytes<R: Rng>(params: &Params, trials: usize, rng: &mut R) -> usize {
    let mut failures = 0;
    let ((a, t), s) = keygen(params, rng);
    for _ in 0..trials {
        let mut input_bytes = Vec::new();
        for _ in 0..(params.n / 8) {
            let byte = rng.gen_range(0..256) as u8;
            input_bytes.push(byte);
        }
        println!("input_bytes.len {}", input_bytes.len());
        let (u, v) = encrypt_bytes(params, &a, &t, &input_bytes, rng);
        let output_bytes = decrypt_bytes(params, &s, &u, &v);
        if input_bytes != output_bytes {
            failures += 1;
            println!("Failure:\ninput:  {input_bytes:?}\noutput: {output_bytes:?}",);
        }
    }
    failures
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let mut rng = thread_rng();
        let params = Params {
            n: 256,
            d: 2,
            q: 65535, // Increased q for more noise tolerance
            eta: 1,
        };
        let trials = 100;
        let failures = test_encryption(&params, trials, &mut rng);
        assert_eq!(failures, 0);
        println!("Failures: {}/{}", failures, trials);
    }

    #[test]
    fn encrypt_decrypt_bytes() {
        let mut rng = thread_rng();
        let params = Params {
            n: 256,
            d: 2,
            q: 65535, // Increased q for more noise tolerance
            eta: 1,
        };
        let trials = 100;
        let failures = test_encrypt_bytes(&params, trials, &mut rng);
        assert_eq!(failures, 0);
        println!("Failures: {}/{}", failures, trials);
    }
}
