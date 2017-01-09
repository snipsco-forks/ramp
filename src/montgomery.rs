
use int::Int;

pub struct ModInt(Int);

#[derive(Debug)]
pub struct Modulus<'a> {
    pub modulus: &'a Int,
    pub modulus_inv0: usize,
    pub limbs: usize,
    pub r: Int,
}

impl<'a> Modulus<'a> {
    #[allow(dead_code)]
    pub fn new(modulus: &'a Int) -> Modulus<'a> {
        use ll::limb::Limb;
        let limbs_count = (modulus.bit_length() as usize + Limb::BITS - 1) / Limb::BITS;
        let r = Int::one() << (limbs_count * Limb::BITS);
        Modulus {
            modulus: modulus,
            modulus_inv0: ::ll::montgomery::single_limb_montgomery_inverse((&r - modulus)
                .limbs()
                .0 as _),
            limbs: limbs_count,
            r: r.clone(),
        }
    }

    fn redc(&self, a: &mut Int) {
        unsafe {
            let mut t = Int::with_capacity(2 * self.limbs as u32);
            ::ll::copy_incr(a.limbs(), t.limbs_uninit(), a.abs_size());
            for i in a.abs_size()..(t.cap as i32) {
                *t.limbs_uninit().offset(i as isize) = ::ll::limb::Limb(0);
            }
            ::ll::montgomery::montgomery_redc(a.limbs_uninit(),
                                              self.limbs as i32,
                                              self.modulus.limbs(),
                                              self.modulus_inv0,
                                              t.limbs_uninit());
            a.size = self.limbs as i32;
        }
    }

    pub fn mul(&self, a: &ModInt, b: &ModInt) -> ModInt {
        unsafe {
            let mut t = Int::with_capacity(2 * self.limbs as u32);
            t.size = t.cap as i32;
            ::ll::mul(t.limbs_uninit(), a.0.limbs(), self.limbs as i32, b.0.limbs(), self.limbs as i32);
            self.redc(&mut t);
            ModInt(t)
        }
    }

    pub fn sqr(&self, a: &ModInt) -> ModInt {
        unsafe {
            let mut t = Int::with_capacity(2 * self.limbs as u32);
            t.size = t.cap as i32;
            ::ll::sqr(t.limbs_uninit(), a.0.limbs(), self.limbs as i32);
            self.redc(&mut t);
            ModInt(t)
        }
    }

    pub fn pow(&self, a: &ModInt, b: &Int) -> ModInt {
        let mut result = self.to_montgomery(&Int::one());
        unsafe {
            ::ll::montgomery::modpow_by_montgomery(result.0.limbs_uninit(),
                                                   self.limbs as i32,
                                                   self.modulus.limbs(),
                                                   a.0.limbs(),
                                                   b.limbs(),
                                                   b.abs_size());
        }
        result
    }

    fn montgomerize(&self, a: &mut Int) {
        Self::pad_to(a, self.limbs);
    }

    fn pad_to(a: &mut Int, s:usize) {
        unsafe {
            a.ensure_capacity(s as u32);
            for i in a.abs_size()..(a.cap as i32) {
                *a.limbs_uninit().offset(i as isize) = ::ll::limb::Limb(0);
            }
            a.size = s as i32;
        }
    }

    #[allow(dead_code)]
    pub fn to_montgomery(&self, a: &Int) -> ModInt {
        let mut it = (a * &self.r) % self.modulus;
        self.montgomerize(&mut it);
        ModInt(it)
    }

    #[allow(dead_code)]
    pub fn to_natural(&self, a: ModInt) -> Int {
        let mut it = a.0;
        it.normalize();
        it %= self.modulus;
        Self::pad_to(&mut it, 2*self.limbs);
        self.redc(&mut it);
        it.normalize();
        it
    }
}

#[cfg(test)]
mod test {
    use ::int::Int;

    #[test]
    fn test_montgomery_redc() {
        let cases = [("1547425065876476735897735405", "193514046488575", "87960930698705")];
        for &(a_bar, m, x_bar) in &cases {
            let mut a_bar = a_bar.parse().unwrap();
            let m = m.parse().unwrap();
            let x_bar: Int = x_bar.parse().unwrap();
            let mg = super::Modulus::new(&m);
            mg.redc(&mut a_bar);
            assert_eq!(a_bar, x_bar);
        }
    }

    #[test]
    fn test_montgomery_cvt() {
        let cases = [("7", "9"), ("6", "4053222090678603523540592804780123937619987201526761")];
        for &(a, m) in &cases {
            let a = a.parse().unwrap();
            let m = m.parse().unwrap();
            let mg = super::Modulus::new(&m);
            assert_eq!(mg.to_natural(mg.to_montgomery(&a)), a);
        }
    }

    #[test]
    fn test_montgomery_mul() {
        let cases = [("2", "13", "207"), ("5", "1", "193514046488575")];
        for &(a, b, m) in &cases {
            let a = a.parse().unwrap();
            let b = b.parse().unwrap();
            let m = m.parse().unwrap();
            let mg = super::Modulus::new(&m);
            let a_bar = mg.to_montgomery(&a);
            let b_bar = mg.to_montgomery(&b);
            let ab_bar = mg.mul(&a_bar, &b_bar);
            let ab = mg.to_natural(ab_bar);
            assert_eq!(ab, a * b % &m);
        }
    }

    // #[test]
    // fn test_montgomery() {
    // let (p, q, n, x) = parse_them();
    // let mut m = super::Montgomery::new(&n);
    // assert_eq!(m.to_natural(&m.to_montgomery(&p)), p);
    // assert_eq!(m.to_natural(&m.to_montgomery(&q)), q);
    // let mut mp = m.to_montgomery(&p);
    // let mq = m.to_montgomery(&q);
    // m.mul(&mut mp, &mq);
    // assert_eq!(m.to_natural(&mp), (&p * &q) % &n);
    // let mut mp = m.to_montgomery(&p);
    // m.sq(&mut mp);
    // assert_eq!(m.to_natural(&mp), (&p * &p) % &n);
    //        assert_eq!(super::modpow_lr_k_ary_with_k_mgmy(&p, &q, &n, 7), x);
    // }
    //

}
