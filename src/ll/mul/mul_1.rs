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
use ll::limb_ptr::{Limbs, LimbsMut};

#[inline(always)]
pub unsafe fn mul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    if cfg!(feature="asm") {
        asm(wp, xp, n, vl)
    } else {
        generic_unroll_4(wp, xp, n, vl)
    }
}

#[allow(dead_code)]
#[inline(always)]
pub unsafe fn generic(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
    let mut cl = Limb(0);
    loop {
        let xl = *xp;
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        *wp = lpl;

        n -= 1;
        if n == 0 { break; }

        wp = wp.offset(1);
        xp = xp.offset(1);
    }

    return cl;
}

#[allow(dead_code)]
#[inline(always)]
pub unsafe fn generic_unroll_4(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
    let mut cl = Limb(0);
    while n / 4 != 0 {
        let xl = *xp;
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;
        *wp = lpl;

        let xl = *xp.offset(1);
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;
        *wp.offset(1) = lpl;

        let xl = *xp.offset(2);
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;
        *wp.offset(2) = lpl;

        let xl = *xp.offset(3);
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;
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
        *wp = lpl;

        n -= 1;
        wp = wp.offset(1);
        xp = xp.offset(1);
    }

    return cl;
}

/**
 * Multiplies the `n` least-significant limbs of `xp` by `vl` storing the `n` least-significant
 * limbs of the product in `{wp, n}`.
 *
 * Returns the highest limb of the product
 */
#[cfg(asm)]
#[inline]
pub unsafe fn asm(mut wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(::ll::same_or_incr(wp, n, xp, n));
    extern "C" {
        fn ramp_mul_1(wp: *mut Limb, xp: *const Limb, n: i32, vl: Limb) -> Limb;
    }

    ramp_mul_1(&mut *wp, &*xp, n, vl)
}

/**
 * Multiplies the `n` least-significant limbs of `xp` by `vl` storing the `n` least-significant
 * limbs of the product in `{wp, n}`.
 *
 * Returns the highest limb of the product
 */
#[cfg(not(asm))]
#[inline]
pub unsafe fn asm(_: LimbsMut, _: Limbs, _: i32, _: Limb) -> Limb {
    unimplemented!()
}

