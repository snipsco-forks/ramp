// Copyright 2015 The Ramp Developers
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

//! Multiprecision Montgomery arithmetics.
//!
//! This module contains optimized modular arithmetic operations under
//! [Montgomery transformation](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication).
//!
//! These optimizations bring significant performance improvements to the
//! modular exponentiation algorithms, or when many operations are performed
//! on the same data within a constant modular field.
//!

use int::Int;

/// A Montgomery modulus.
///
/// This structure holds precomputed values that optimized subsequent
/// computation on MtgyInt.
///
/// It can only be constructed for positive odd Modulus.
///
/// # Examples
///
/// Starting with 17 as a modulus, prepare the modulus helpers.
///
/// ```rust
/// use ramp::int::Int;
/// use ramp::int::mtgy::*;
///
/// let m:Int = 17.into();
/// let modulus = MtgyModulus::new(&m);
/// ```
///
/// Convert between Montgomery and natural space:
///
/// ```rust
/// # use ramp::int::Int;
/// # use ramp::int::mtgy::*;
/// # let m:Int = 17.into();
/// # let modulus = MtgyModulus::new(&m);
/// let a:Int = 5.into();
/// let a_bar = modulus.to_mtgy(&a);
/// let a_back = modulus.to_int(&a_bar);
/// assert_eq!(a, a_back);
/// ```
///
/// Perform a modular multiplication in Montgomery space:
///
/// ```rust
/// # use ramp::int::Int;
/// # use ramp::int::mtgy::*;
/// # let m:Int = 17.into();
/// # let modulus = MtgyModulus::new(&m);
/// let a:Int = 5.into();
/// let a_bar = modulus.to_mtgy(&a);
/// let b:Int = 7.into();
/// let b_bar = modulus.to_mtgy(&b);
/// let ab_bar = modulus.mul(&a_bar, &b_bar);
/// let ab = modulus.to_int(&ab_bar);
/// assert_eq!(ab, a*b % &m);
/// ```
///
/// Perform a modular exponentiation in Montgomery space. Note that while the
/// basis is in Montgomery form, the exponent is in natural space.
///
/// ```rust
/// # use ramp::int::Int;
/// # use ramp::int::mtgy::*;
/// # let m:Int = 17.into();
/// # let modulus = MtgyModulus::new(&m);
/// let a:Int = 5.into();
/// let a_bar = modulus.to_mtgy(&a);
/// let a_pow_7_bar = modulus.pow(&a_bar, &Int::from(7));
/// let a_pow_7 = modulus.to_int(&a_pow_7_bar);
/// assert_eq!(a_pow_7, a.pow(7) % &m);
/// ```
///
#[derive(Debug)]
pub struct MtgyModulus<'a> {
    modulus: &'a Int,
    modulus_inv0: ::ll::limb::Limb,
    limbs: usize,
    r: Int,
}

/// An integer in Montgomery form.
///
/// The Montgomery form is valid for one and only one MtgyModulus. It's the
/// user responsibility to maintain this consistency (aka, don't mix up
/// MtgyInt from different MtgyModulus).
pub struct MtgyInt(Int);

impl<'a> MtgyModulus<'a> {
    /// Builds a pre-optimized MtgyModulus to perform.
    ///
    /// # Panic
    ///
    /// For the Montgomery form to exists, the modulus has to be odd (and positive).
    /// The constructor will panic otherwise.
    #[allow(dead_code)]
    pub fn new(modulus: &'a Int) -> MtgyModulus<'a> {
        assert!(!modulus.is_even(), "Montgomery modulus must be odd");
        assert_eq!(modulus.sign(), 1, "Montgomery modulus must be positive");
        use ll::limb::Limb;
        let limbs_count = (modulus.bit_length() as usize + Limb::BITS - 1) / Limb::BITS;
        let r = Int::one() << (limbs_count * Limb::BITS);
        MtgyModulus {
            modulus: modulus,
            modulus_inv0: ::ll::mtgy::inv1(*(&r - modulus).limbs()),
            limbs: limbs_count,
            r: r.clone(),
        }
    }

    fn redc(&self, a: &mut Int) {
        unsafe {
            assert_eq!(a.abs_size(), 2*self.limbs as i32);
            let mut t = Int::with_capacity(2 * self.limbs as u32);
            ::ll::copy_incr(a.limbs(), t.limbs_uninit(), a.abs_size());
            for i in a.abs_size()..(t.cap as i32) {
                *t.limbs_uninit().offset(i as isize) = ::ll::limb::Limb(0);
            }
            ::ll::mtgy::redc(a.limbs_uninit(),
                             self.limbs as i32,
                             self.modulus.limbs(),
                             self.modulus_inv0,
                             t.limbs_uninit());
            a.size = self.limbs as i32;
        }
    }

    /// Multiply two integers under Montgomery form.
    ///
    /// # Panic
    ///
    /// Panics if the two integers are not of the expected size (it is
    /// only likely to happen in case of a mixup of two MtgyModulus).
    pub fn mul(&self, a: &MtgyInt, b: &MtgyInt) -> MtgyInt {
        unsafe {
            assert_eq!(a.0.abs_size(), self.limbs as i32);
            assert_eq!(b.0.abs_size(), self.limbs as i32);
            let mut t = Int::with_capacity(2 * self.limbs as u32);
            t.size = t.cap as i32;
            ::ll::mul(t.limbs_uninit(),
                      a.0.limbs(),
                      self.limbs as i32,
                      b.0.limbs(),
                      self.limbs as i32);
            self.redc(&mut t);
            MtgyInt(t)
        }
    }

    /// Square an integer in Montgomery form.
    ///
    /// # Panic
    ///
    /// Panics if the integer is not of the expected size (it is
    /// only likely to happen in case of a mixup of two MtgyModulus).
    pub fn sqr(&self, a: &MtgyInt) -> MtgyInt {
        unsafe {
            assert_eq!(a.0.abs_size(), self.limbs as i32);
            let mut t = Int::with_capacity(2 * self.limbs as u32);
            t.size = t.cap as i32;
            ::ll::sqr(t.limbs_uninit(), a.0.limbs(), self.limbs as i32);
            self.redc(&mut t);
            MtgyInt(t)
        }
    }

    /// Compute a modular exponentiation under Montgomery form.
    ///
    /// Note that `basis` is expected in Montgomery form, while `exponent` 
    /// is a natural int.
    ///
    /// # Panic
    ///
    /// * Panics if the basis integer is not of the expected size (it is
    /// only likely to happen in case of a mixup of two MtgyModulus).
    /// * Panics if exponent is negative.
    pub fn pow(&self, basis: &MtgyInt, exponent: &Int) -> MtgyInt {
        let mut result = self.to_mtgy(&Int::one());
        unsafe {
            assert_eq!(basis.0.abs_size(), self.limbs as i32);
            assert!(exponent.sign() >= 0);
            ::ll::mtgy::modpow(result.0.limbs_uninit(),
                               self.limbs as i32,
                               self.modulus.limbs(),
                               self.modulus_inv0,
                               basis.0.limbs(),
                               exponent.limbs(),
                               exponent.abs_size());
        }
        result
    }

    fn montgomerize(&self, a: &mut Int) {
        Self::pad_to(a, self.limbs);
    }

    fn pad_to(a: &mut Int, s: usize) {
        unsafe {
            a.ensure_capacity(s as u32);
            for i in a.abs_size()..(a.cap as i32) {
                *a.limbs_uninit().offset(i as isize) = ::ll::limb::Limb(0);
            }
            a.size = s as i32;
        }
    }

    /// Convert an int to its Montgomery form.
    #[allow(dead_code)]
    pub fn to_mtgy(&self, a: &Int) -> MtgyInt {
        let mut it = (a * &self.r) % self.modulus;
        self.montgomerize(&mut it);
        MtgyInt(it)
    }

    /// Convert a Montgomery int back to Int.
    /// # Panic
    ///
    /// * Panics if the integer is not of the expected size (it is
    /// only likely to happen in case of a mixup of two MtgyModulus).
    #[allow(dead_code)]
    pub fn to_int(&self, a: &MtgyInt) -> Int {
        assert_eq!(a.0.abs_size(), self.limbs as i32);
        let mut it = unsafe {
            let mut it = Int::with_capacity(2 * self.limbs as u32);
            ::ll::copy_incr(a.0.limbs(), it.limbs_uninit(), self.limbs as i32);
            it.size = self.limbs as i32;
            it.normalize();
            it
        };
        it %= self.modulus;
        Self::pad_to(&mut it, 2 * self.limbs);
        self.redc(&mut it);
        it.normalize();
        it
    }
}

#[cfg(test)]
mod test {
    use ::int::Int;

    #[test]
    fn redc() {
        let cases = [("1547425065876476735897735405", "193514046488575", "87960930698705")];
        for &(a_bar, m, x_bar) in &cases {
            let mut a_bar = a_bar.parse().unwrap();
            let m = m.parse().unwrap();
            let x_bar: Int = x_bar.parse().unwrap();
            let mg = super::MtgyModulus::new(&m);
            mg.redc(&mut a_bar);
            assert_eq!(a_bar, x_bar);
        }
    }

    #[test]
    fn cvt() {
        let cases = [("1", "1009"),
                     ("15", "1009"),
                     ("9330786055998253486590", "4349330786055998253486590232462401"),
                     ("7", "9"),
                     ("6", "4053222090678603523540592804780123937619987201526761")];
        for &(a, m) in &cases {
            let a = a.parse().unwrap();
            let m = m.parse().unwrap();
            let mg = super::MtgyModulus::new(&m);
            assert_eq!(mg.to_int(&mg.to_mtgy(&a)), a);
        }
    }

    #[test]
    fn mul() {
        let cases = [
            ("1", "2", "13", "2"),
            ("1", "1", "13", "1"),
            ("7", "7", "13", "10"),
            ("2", "13", "207", "26"),
            ("1", "1", "1009", "1"),
            ("2", "10", "1009", "20"),
            ("5", "1", "193514046488575", "5"),
            ("15", "1", "4349330786055998253486590232462401", "15"),
            ("15", "10", "1475703270992002140168997557525132617116077748043980354291003276386587324053694848174953095546817655706234979251318204003655882580688895", "150"),
            ("148677972634832330983979593310074301486537017973460461278300587514468301043894574906886127642530475786889672304776052879927627556769456140664043088700743909632312483413393134504352834240399191134336344285483935856491230340093391784574980688823380828143810804684752914935441384845195613674104960646037368551517",
            "158741574437007245654463598139927898730476924736461654463975966787719309357536545869203069369466212089132653564188443272208127277664424448947476335413293018778018615899291704693105620242763173357203898195318179150836424196645745308205164116144020613415407736216097185962171301808761138424668335445923774195463",
            "446397596678771930935753654586920306936946621208913265356418844327220812727766442444894747633541329301877801861589929170469310562024276317335720389819531817915083642419664574530820516411614402061341540773621609718596217130180876113842466833544592377419546315874157443700724565446359813992789873047692473646165446397596678771930935753654586920306936946621208913265356418844327220812727766442444894747633541329301877801861589929170469310562045923774195463",
            "15733033542428556326610775226428250291950090984377467644096837926072\
            98553857572965450727431838091748906310425930542328045644280094594289\
            52380420588404540083723320848855612172087517363909606183916778041064\
            11997952939978862543172484483575568826983703005515400230343351224994\
            85403291437917132468481025327704901371719125205664144192914895118949\
            25716605685210349843822514310138216212323303683754146084454361295646\
            557462263542138176646203699553393662651092450")
        ];
        for &(a, b, m, x) in &cases {
            let a = a.parse().unwrap();
            let b = b.parse().unwrap();
            let m = m.parse().unwrap();
            let x:Int = x.parse().unwrap();
            let mg = super::MtgyModulus::new(&m);
            let a_bar = mg.to_mtgy(&a);
            let b_bar = mg.to_mtgy(&b);
            let ab_bar = mg.mul(&a_bar, &b_bar);
            let ab = mg.to_int(&ab_bar);
            assert_eq!(ab, x);
        }
    }

}
