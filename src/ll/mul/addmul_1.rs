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

#![allow(improper_ctypes)]

use ll::limb::Limb;
use ll::same_or_separate;

use ll::limb_ptr::{Limbs, LimbsMut};

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and adds them to the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result.
 */
#[inline]
pub fn addmul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    if cfg!(feature="asm") {
        unsafe { asm(wp, xp, n, vl) }
    } else {
       unsafe {  generic_unroll_4(wp, xp, n, vl) }
    }
}

#[allow(dead_code)]
#[inline]
pub unsafe fn generic(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_separate(wp, n, xp, n));

    let mut cl = Limb(0);
    loop {
        let xl = *xp;
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp).add_overflow(lpl);
        cl = cl + carry;

        *wp = lpl;

        n -= 1;
        if n == 0 { break; }

        wp = wp.offset(1);
        xp = xp.offset(1);
    }

    return cl;
}

#[allow(dead_code)]
#[inline]
pub unsafe fn generic_unroll_4(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_separate(wp, n, xp, n));

    let mut cl = Limb(0);
    while n / 4 != 0 {
        let xl = *xp;
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp).add_overflow(lpl);
        cl = cl + carry;
        *wp = lpl;

        let xl = *xp.offset(1);
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp.offset(1)).add_overflow(lpl);
        cl = cl + carry;
        *wp.offset(1) = lpl;

        let xl = *xp.offset(2);
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp.offset(2)).add_overflow(lpl);
        cl = cl + carry;
        *wp.offset(2) = lpl;

        let xl = *xp.offset(3);
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp.offset(3)).add_overflow(lpl);
        cl = cl + carry;
        *wp.offset(3) = lpl;

        n -= 4;
        wp = wp.offset(4);
        xp = xp.offset(4);
    }

    while n != 0 {
        let xl = *xp;
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp).add_overflow(lpl);
        cl = cl + carry;

        *wp = lpl;

        n -= 1;
        wp = wp.offset(1);
        xp = xp.offset(1);
    }

    return cl;
}

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and adds them to the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result.
 */
#[cfg(asm)]
#[inline]
pub unsafe fn asm(mut wp: LimbsMut, xp:  Limbs, n: i32, vl: Limb) -> Limb {
    extern "C" {
        fn ramp_addmul_1(wp: *mut Limb, xp: *const Limb, n: i32, vl: Limb) -> Limb;
    }

    ramp_addmul_1(&mut *wp, &*xp, n, vl)
}

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and adds them to the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result.
 */
#[cfg(not(asm))]
#[inline]
pub unsafe fn asm(_: LimbsMut, _: Limbs, _: i32, _: Limb) -> Limb {
    unimplemented!()
}

#[cfg(test)]
mod test {

    macro_rules! t {
        ($func:ident) => {
            #[test]
            fn $func() {
                use ll::limb::Limb;
                use ll::limb_ptr::{Limbs, LimbsMut};
                let half_limb = 1usize << (Limb::BITS-1);
                unsafe {
                    for &(w, x, l, exp_w, exp_c) in &[
                        (&[0usize, 0] as &[usize], &[0usize, 0] as &[usize], 0usize, &[0usize, 0] as &[usize], 0),
                        (&[0, 0], &[1, 0], 1, &[1, 0], 0),
                        (&[0, 0], &[1, 1], 1, &[1, 1], 0),
                        (&[!0, !0], &[1, 0], 1, &[0, 0], 1),
                        (&[0, 0, 0, 0], &[1, 0], 1, &[1, 0, 0, 0], 0),
                        (&[0, 0, 0, 0, 0], &[1, 0], 1, &[1, 0, 0, 0, 0], 0),
                        (&[0, 0, 0, 0, 0], &[1, 2, 3, 4, 5], 1, &[1, 2, 3, 4, 5], 0),
                        (&[0, 0, 0, 0, 0], &[1, 2, 3, 4, 5], 2, &[2, 4, 6, 8, 10], 0),
                        (&[1, 2, 3, 4, 5], &[0, 0, 0, 0, 0], 1, &[1, 2, 3, 4, 5], 0),
                        (&[0, 0, 0, 0, 0], &[1, 2, 3, 4, 5], 2, &[2, 4, 6, 8, 10], 0),
                        (&[6, 7, 8, 9,10], &[1, 2, 3, 4, 5], 3, &[9, 13, 17, 21, 25], 0),
                        (&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1, 
                            &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 0),
                        (&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 1, 
                            &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 0),
                        (&[!0, !0, !0, !0, !0], &[1, 0, 0, 0, 0], 1, &[0, 0, 0, 0, 0], 1),
                        (&[0, 0, 0, 0, 0], &[half_limb, half_limb, half_limb, half_limb, half_limb], 2,
                            &[0, 1, 1, 1, 1], 1),
                    ] {
                        let w_vec = w.to_vec();
                        let x_vec = x.to_vec();
                        let x_limbs = Limbs::new(x_vec.as_ptr() as _, 0, x.len() as i32);
                        let w_limbs = LimbsMut::new(w_vec.as_ptr() as _, 0, w_vec.len() as i32);
                        let Limb(c) = super::$func(w_limbs, x_limbs, x.len() as _, Limb(l as _));
                        assert_eq!(exp_w, &*w_vec,
                                   "wrong result testing {:?}+{:?}*{:?}={:?}", w, x, l, w_vec);
                        assert_eq!(exp_c, c,
                                   "wrong carry testing {:?}+{:?}*{:?}={:?}", w, x, l, w_vec);
                    }
                }
            }
        }
    }

    t!(generic);
    t!(generic_unroll_4);
}
