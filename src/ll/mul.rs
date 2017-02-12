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

use std::cmp::Ordering;

use ll;
use ll::limb::Limb;
use super::{overlap, same_or_separate, same_or_incr};
use mem;

use ll::limb_ptr::{Limbs, LimbsMut};

const TOOM22_THRESHOLD : i32 = 20;

#[allow(dead_code)]
#[inline]
unsafe fn mul_1_generic(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
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

/**
 * Multiplies the `n` least-significant limbs of `xp` by `vl` storing the `n` least-significant
 * limbs of the product in `{wp, n}`.
 *
 * Returns the highest limb of the product
 */
#[inline]
//#[cfg(not(target_arch="x86_64"))]
pub unsafe fn mul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_incr(wp, n, xp, n));

    mul_1_generic(wp, xp, n, vl)
}

/**
 * Multiplies the `n` least-significant limbs of `xp` by `vl` storing the `n` least-significant
 * limbs of the product in `{wp, n}`.
 *
 * Returns the highest limb of the product
 */
#[inline]
#[cfg(target_arch="x86_64")]
#[allow(unused_assignments)]
pub unsafe fn _mul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_incr(wp, n, xp, n));
    let mut r:usize = 0;
    let mut n:i64 = n as _;
    let mut w:*mut _ = &mut *wp.offset(0);
    let mut x:*const _ = &*xp.offset(0);
    while n % 4 != 0 {
        asm!("
        movq ($2), %rax
        movq %rdx, %r8
        mulq $8
        addq %r8, %rax
        adcq $$0, %rdx
        movq %rax, ($1)
        add $$8, $2
        add $$8, $1
        sub $$1, $3
        "
        : "=&{rdx}"(r), "=&r"(w), "=&r"(x), "=&r"(n)
        : "0"(r), "1"(w), "2"(x), "3"(n), "r"(vl.0)
        : "r8", "rax", "memory", "cc");
    }
    if n != 0 {
        asm!("
        lea ($1,$3,8), $1
        lea ($2,$3,8), $2
        neg $3

        .align 4
        1: 
        mov ($2,$3,8), %rax
        mul $8
        add %rax, %r8
        adc $$0, %rdx
        mov 8($2,$3,8), %rax
        mov %r8, ($1,$3,8)
        mov %rdx, %r9

        mul $8
        add %rax, %r9
        adc $$0, %rdx
        mov 16($2,$3,8), %rax
        mov %r9, 8($1,$3,8)
        mov %rdx, %r10

        mul $8
        add %rax, %r10
        adc $$0, %rdx
        mov 24($2,$3,8), %rax
        mov %r10, 16($1,$3,8)
        mov %rdx, %r11

        mul $8
        add %rax, %r11
        adc $$0, %rdx
        mov %rdx, %r8
        mov %r11, 24($1,$3,8)

        add $$4, $3
        js 1b
        "
        : "=&{r8}"(r), "=&r"(w), "=&r"(x), "=&r"(n)
        : "0"(r), "1"(w), "2"(x), "3"(n), "r"(vl.0)
        : "r9", "r10", "r11", "rax", "rdx", "memory", "cc");
    }
    Limb(r as _)
}

#[inline]
#[allow(dead_code)]
unsafe fn addmul_1_generic(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
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

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and adds them to the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result.
 */
#[inline]
//#[cfg(not(target_arch="x86_64"))]
pub unsafe fn addmul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    addmul_1_generic(wp, xp, n, vl)
}

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and adds them to the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result.
 */
#[inline]
#[cfg(target_arch="x86_64")]
#[allow(unused_assignments)]
pub unsafe fn _addmul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_incr(wp, n, xp, n));
    let mut r:usize = 0;
    let mut n:i64 = n as _;
    let mut w:*mut _ = &mut *wp.offset(0);
    let mut x:*const _ = &*xp.offset(0);
    while n % 4 != 0 {
        asm!("
        mov ($2), %rax
        mov %rdx, %r8
        mul $8
        add %r8, %rax
        adc $$0, %rdx
        add %rax, ($1)
        adc $$0, %rdx
        add $$8, $2
        add $$8, $1
        sub $$1, $3
        "
        : "=&{rdx}"(r), "=&r"(w), "=&r"(x), "=&r"(n)
        : "0"(r), "1"(w), "2"(x), "3"(n), "r"(vl.0)
        : "r8", "rax", "memory", "cc");
    }
    if n != 0 {
        asm!("
        lea ($1,$3,8), $1
        lea ($2,$3,8), $2
        neg $3

        .align 4
        1: 
        mov ($2,$3,8), %rax
        mul $8
        add %rax, %r8
        adc $$0, %rdx
        mov 8($2,$3,8), %rax
        add %r8, ($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r9

        mul $8
        add %rax, %r9
        adc $$0, %rdx
        mov 16($2,$3,8), %rax
        add %r9, 8($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r10

        mul $8
        add %rax, %r10
        adc $$0, %rdx
        mov 24($2,$3,8), %rax
        add %r10, 16($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r11

        mul $8
        add %rax, %r11
        adc $$0, %rdx
        add %r11, 24($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r8

        add $$4, $3
        js 1b
        "
        : "=&{r8}"(r), "=&r"(w), "=&r"(x), "=&r"(n)
        : "0"(r), "1"(w), "2"(x), "3"(n), "r"(vl.0)
        : "r9", "r10", "r11", "rax", "rdx", "memory", "cc");
    }
    Limb(r as _)
}
/*
#[inline(always)]
#[allow(dead_code)]
pub unsafe fn _addmul_2(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl1: Limb, vl2:Limb) -> (Limb, Limb) {
    debug_assert!(n > 0);
    debug_assert!(same_or_separate(wp, n, xp, n));

    let mut carry_1 = Limb(0);
    let mut carry_2 = Limb(0);
    loop {
        let (hi_l2, lo_l2) = (*xp).mul_hilo(vl2);
        let (hi_l1, lo_l1) = (*xp).mul_hilo(vl1);

        let (s, cl1) = (*wp).add_overflow(lo_l1);
        let (s, cl2) = s.add_overflow(carry_1);
        *wp = s;

        let (s, cl1) = carry_2.add_overflow(Limb(cl1 as _) + Limb(cl2 as _));
        let (s, cl2) = s.add_overflow(lo_l2);
        let (s, cl3) = s.add_overflow(hi_l1);
        carry_1 = s;

        carry_2 = hi_l2 + (Limb(cl1 as _) + Limb(cl2 as _) + Limb(cl3 as _));

        n -= 1;
        if n == 0 { break; }
        xp = xp.offset(1);
        wp = wp.offset(1);
    }
    (carry_1, carry_2)
}

#[inline(always)]
#[allow(dead_code)]
#[allow(unused_assignments)]
pub unsafe fn addmul_2(wp: LimbsMut, xp: Limbs, n: i32, vl1: Limb, vl2:Limb) -> (Limb, Limb) {
    debug_assert!(n > 0);
    debug_assert!(same_or_separate(wp, n, xp, n));
    let mut n:i64 = n as _;
    let mut w:*mut _ = &mut *wp.offset(0);
    let mut x:*const _ = &*xp.offset(0);

    let mut carry_1 = 0u64;
    let mut carry_2 = 0u64;
    asm!("
    mov ($3), %rbx
    lea ($2,$4,8), $2       // $2 is wp
    lea ($3,$4,8), $3       // $3 is xp
    neg $4                  // $4 is n

    .align 4
    1: 
    xor %r8, %r8
    mov 8($3,$4,8), %rbx
    mov %rbx, %rax
    mul $10                 // $10 is vl1
    add %rax, $0
    adc %rdx, $1
    adc $$0, %r8            // temps: %r8:$1:$0
    add $0, ($2,$4,8)       // write wi
    adc $$0, $1
    adc $$0, %r8
    mov $1, $0
    mov %r8, $1
    mov %rbx, %rax
    mov 8($3,$4,8), %rbx
    mul $11                 // $11 is vl2
    add %rax, $0
    adc %rdx, $1

    add $$1, $4
    js 1b
    "
    : "=&r"(carry_1), "=&r"(carry_2), "=&r"(w), "=&r"(x), "=&r"(n)
    : "0"(carry_1), "1"(carry_2), "2"(w), "3"(x), "4"(n), "r"(vl1.0), "r"(vl2.0)
    : "r8", "rax", "rbx", "rdx", "memory", "cc");
    (Limb(carry_1), Limb(carry_2))
}
*/

#[inline]
#[allow(dead_code)]
unsafe fn submul_1_generic(mut wp: LimbsMut, mut xp: Limbs, mut n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_separate(wp, n, xp, n));

    let mut cl = Limb(0);
    loop {
        let xl = *xp;
        let (hpl, lpl) = xl.mul_hilo(vl);
        let (lpl, carry) = lpl.add_overflow(cl);
        cl = hpl + carry;

        let (lpl, carry) = (*wp).sub_overflow(lpl);
        cl = cl + carry;

        *wp = lpl;

        n -= 1;
        if n == 0 { break; }

        wp = wp.offset(1);
        xp = xp.offset(1);
    }

    return cl;
}

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and subtracts them from the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result, adjust for borrow.
 */
#[cfg(not(asm))]
#[inline]
#[cfg(not(target_arch="x86_64"))]
pub unsafe fn submul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    submul_1_generic(wp, xp, n, vl)
}

/**
 * Multiplies the `n` least-signficiant digits of `xp` by `vl` and subtracts them from the `n`
 * least-significant digits of `wp`. Returns the highest limb of the result, adjust for borrow.
 */
#[inline]
#[cfg(target_arch="x86_64")]
#[allow(unused_assignments)]
pub unsafe fn submul_1(wp: LimbsMut, xp: Limbs, n: i32, vl: Limb) -> Limb {
    debug_assert!(n > 0);
    debug_assert!(same_or_incr(wp, n, xp, n));
    let mut r:usize = 0;
    let mut n:i64 = n as _;
    let mut w:*mut _ = &mut *wp.offset(0);
    let mut x:*const _ = &*xp.offset(0);
    while n % 4 != 0 {
        asm!("
        mov ($2), %rax
        mov %rdx, %r8
        mul $8
        add %r8, %rax
        adc $$0, %rdx
        sub %rax, ($1)
        adc $$0, %rdx
        add $$8, $2
        add $$8, $1
        sub $$1, $3
        "
        : "=&{rdx}"(r), "=&r"(w), "=&r"(x), "=&r"(n)
        : "0"(r), "1"(w), "2"(x), "3"(n), "r"(vl.0)
        : "r8", "rax", "memory", "cc");
    }
    if n != 0 {
        asm!("
        lea ($1,$3,8), $1
        lea ($2,$3,8), $2
        neg $3

        .align 4
        1: 
        mov ($2,$3,8), %rax
        mul $8
        add %rax, %r8
        adc $$0, %rdx
        mov 8($2,$3,8), %rax
        sub %r8, ($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r9

        mul $8
        add %rax, %r9
        adc $$0, %rdx
        mov 16($2,$3,8), %rax
        sub %r9, 8($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r10

        mul $8
        add %rax, %r10
        adc $$0, %rdx
        mov 24($2,$3,8), %rax
        sub %r10, 16($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r11

        mul $8
        add %rax, %r11
        adc $$0, %rdx
        sub %r11, 24($1,$3,8)
        adc $$0, %rdx
        mov %rdx, %r8

        add $$4, $3
        js 1b
        "
        : "=&{r8}"(r), "=&r"(w), "=&r"(x), "=&r"(n)
        : "0"(r), "1"(w), "2"(x), "3"(n), "r"(vl.0)
        : "r9", "r10", "r11", "rax", "rdx", "memory", "cc");
    }
    Limb(r as _)
}

/**
 * Multiplies `{xp, xs}` by `{yp, ys}`, storing the result to `{wp, xs + ys}`.
 *
 * `{wp, xs + ys}` must be disjoint from both inputs.
 */
pub unsafe fn mul(wp: LimbsMut, xp: Limbs, xs: i32, yp: Limbs, ys: i32) {
    debug_assert!(xs >= ys);
    debug_assert!(ys > 0);
    debug_assert!(!overlap(wp, xs + ys, xp, xs));
    debug_assert!(!overlap(wp, xs + ys, yp, ys));

    // TODO: Pick between algorithms based on input sizes
    if ys <= TOOM22_THRESHOLD {
        mul_basecase(wp, xp, xs, yp, ys);
    } else {
        let mut tmp = mem::TmpAllocator::new();
        let scratch = tmp.allocate((xs * 2) as usize);

        // Can't use xs >= (ys * 2) because if xs is odd, some other invariants
        // in toom22 don't hold
        if (xs * 2) >= (ys * 3) {
            mul_unbalanced(wp, xp, xs, yp, ys, scratch);
        } else {
            mul_toom22(wp, xp, xs, yp, ys, scratch);
        }
    }
}

#[inline(always)]
unsafe fn mul_basecase(mut wp: LimbsMut, xp: Limbs, xs: i32, mut yp: Limbs, mut ys: i32) {
    *wp.offset(xs as isize) = ll::mul_1(wp, xp, xs, *yp);
    wp = wp.offset(1);
    yp = yp.offset(1);
    ys -= 1;

    while ys > 0 {
        *wp.offset(xs as isize) = ll::addmul_1(wp, xp, xs, *yp);
        wp = wp.offset(1);
        yp = yp.offset(1);
        ys -= 1;
    }
}

/*
#[inline(always)]
unsafe fn mul_basecase_addmul_2(mut wp: LimbsMut, xp: Limbs, xs: i32, mut yp: Limbs, mut ys: i32) {
    *wp.offset(xs as isize) = ll::mul_1(wp, xp, xs, *yp);
    wp = wp.offset(1);
    yp = yp.offset(1);
    ys -= 1;

    while ys > 1 {
        let (c1, c2) = addmul_2(wp, xp, xs, *yp, *yp.offset(1));
        *wp.offset(xs as isize) = c1;
        *wp.offset(xs as isize + 1) = c2;

        wp = wp.offset(2);
        yp = yp.offset(2);
        ys -= 2;
    }
    if ys == 1 {
        *wp.offset(xs as isize) = ll::addmul_1(wp, xp, xs, *yp);
    }
}
*/
// Helper fn
#[inline(always)]
pub unsafe fn mul_rec(wp: LimbsMut,
           xp: Limbs, xs: i32,
           yp: Limbs, ys: i32,
           scratch: LimbsMut) {
    if ys < TOOM22_THRESHOLD {
        mul_basecase(wp, xp, xs, yp, ys);
    } else if (xs * 2) >= (ys*3) {
        mul_unbalanced(wp, xp, xs, yp, ys, scratch);
    } else {
        mul_toom22(wp, xp, xs, yp, ys, scratch);
    }
}

unsafe fn mul_toom22(wp: LimbsMut,
                     xp: Limbs, xs: i32,
                     yp: Limbs, ys: i32,
                     scratch: LimbsMut) {
    // Split x into x1, x0 where x = x1*(B^n) + x0
    // Split y into y1, y0 where y = y1*(B^n) + y0
    //
    // Which means the following holds:
    //
    //    x*y = (B^2n + B^n)*z2 - (B^n)*z1 + (B^n + 1)*z0
    //        = (B^2n)*z2 + (B^n)*(z2 + z0 - z1) + z0
    //
    // Where:
    //   z0 = x0*y0
    //   z1 = (x1-x0)*(y1-y0)
    //   z2 = x1*y1
    //
    // z1 is split further into:
    //
    //  zx1 = x1-x0
    //  zy1 = y1-y0
    //
    // So z1 = zx1*zy1

    debug_assert!(xs >= ys && xs < ys*2,
                  "assertion failed: `xs >= ys && xs < ys*2` xs: {}, ys: {}", xs, ys);

    let xh = xs >> 1; // Number of high limbs in x
    let nl = xs - xh; // Number of low limbs
    let yh = ys - nl; // Number of high limbs in y

    debug_assert!(0 < xh && xh <= nl);
    debug_assert!(0 < yh && yh <= xh,
                  "assertion failed: 0 < yh && yh <= xh, xs: {}, ys: {}, xh: {}, yh: {}",
                  xs, ys, xh, yh);

    let x0 = xp; // nl limbs
    let y0 = yp; // nl limbs

    let x1 = xp.offset(nl as isize); // xh limbs
    let y1 = yp.offset(nl as isize); // yh limbs

    let zx1 = wp; // nl limbs
    let zy1 = wp.offset(nl as isize); // nl limbs
    let mut z1_neg = false; // Keep track of whether the real z1 is negative

    // Calculate zx1
    if nl == xh {
        if ll::cmp(x0, x1, nl) == Ordering::Less {
            ll::sub_n(zx1, x1, x0, nl);
            z1_neg = true;
        } else {
            ll::sub_n(zx1, x0, x1, nl);
        }
    } else { // nl > xh
        if ll::is_zero(x0.offset(xh as isize), nl-xh) && ll::cmp(x0, x1, xh) == Ordering::Less {
            ll::sub_n(zx1, x1, x0, xh);
            ll::zero(zx1.offset(xh as isize), nl-xh); // Zero the extra limbs
            z1_neg = true;
        } else {
            ll::sub(zx1, x0, nl, x1, xh);
        }
    }

    // Calculate zy1
    if nl == yh {
        if ll::cmp(y0, y1, nl) == Ordering::Less {
            ll::sub_n(zy1, y1, y0, nl);
            z1_neg = !z1_neg;
        } else {
            ll::sub_n(zy1, y0, y1, nl);
        }
    } else { // nl > yh
        if ll::is_zero(y0.offset(yh as isize), nl-yh) && ll::cmp(y0, y1, yh) == Ordering::Less {
            ll::sub_n(zy1, y1, y0, yh);
            ll::zero(zy1.offset(yh as isize), nl-yh); // Zero the extra limbs
            z1_neg = !z1_neg;
        } else {
            ll::sub(zy1, y0, nl, y1, yh);
        }
    }

    let z0 = wp;
    let z1 = scratch;
    let z2 = wp.offset((nl * 2) as isize);
    let scratch_out = scratch.offset((nl * 2) as isize);

    // Calculate z1 - 2*nl limbs
    mul_rec(z1, zx1.as_const(), nl, zy1.as_const(), nl, scratch_out);

    // Calculate z0 - 2*nl limbs
    mul_rec(z0, x0, nl, y0, nl, scratch_out);

    // Calculate z2 - xh+yh limbs
    mul_rec(z2, x1, xh, y1, yh, scratch_out);

    // Now {wp, 2*nl} = z0 and {wp + 2*nl, xh+yh} = z2

    // {wp + nl, 2*nl} += z0 + z2 - z1
    //                 += {wp, 2*nl}
    //                  + {wp + 2*nl, xh+yh}
    //                  - z1
    //
    // So with {wp, xs+ys}:
    //
    // 0        nl      2*nl        xs+ys
    // +--------+--------+--------+---+
    // |       z0        |     z2     |
    // +--------+--------+--------+---+
    //   +      |        z0       |
    //          +--------+---+----+
    //   +      |       z2   |
    //          +------------+
    //
    // {wp + nl, nl} = HIGH(z0) + LOW(z0) + LOW(z2)
    // {wp + 2*nl, nl} = HIGH(z0) + HIGH(z2) + LOW(z2) + carry

    // LOW(z2) = HIGH(z0) + LOW(z2)
    let cy = ll::add_n(wp.offset((2*nl) as isize),
                       z2.as_const(), z0.offset(nl as isize).as_const(),
                       nl);

    // LOW(z0) + LOW(z2)
    let cy2 = cy + ll::add_n(wp.offset(nl as isize),
                             z0.as_const(), z2.as_const(),
                             nl);

    // LOW(z2) + HIGH(z2)
    let mut cy = cy + ll::add(wp.offset((2*nl) as isize),
                              z2.as_const(), nl,
                              z2.offset(nl as isize).as_const(), yh+xh-nl);

    // Add or subtract `z1` depending on the sign of the real result
    // (we calculate such that it's always positive, easier this way)
    if z1_neg {
        cy = cy + ll::add_n(wp.offset(nl as isize),
                            wp.offset(nl as isize).as_const(), z1.as_const(),
                            2*nl);
    } else {
        cy = cy - ll::sub_n(wp.offset(nl as isize),
                            wp.offset(nl as isize).as_const(), z1.as_const(),
                            2*nl);
    }

    // Apply the carries, has to be done last.
    ll::incr(wp.offset((nl * 2) as isize), cy2);
    ll::incr(wp.offset((nl * 3) as isize), cy);
}

/**
 * Handles multiplication when xs is much bigger than ys.
 *
 * Works basically the same way `mul_1` does, except with `ys` limbs
 * instead of a single limb.
 */
unsafe fn mul_unbalanced(mut wp: LimbsMut,
                         mut xp: Limbs, mut xs: i32,
                         yp: Limbs, ys: i32,
                         scratch: LimbsMut) {
    debug_assert!(xs  > ys);

    mul_toom22(wp, xp, ys, yp, ys, scratch);

    xs -= ys;
    xp = xp.offset(ys as isize);
    wp = wp.offset(ys as isize);

    // Temporary storage for the output of the multiplication
    // in the loop, the loop only needs ys*2 limbs, but the last
    // multiplication needs slightly more than that, but no more
    // than ys*3
    let mut tmp = mem::TmpAllocator::new();
    let w_tmp = tmp.allocate((ys * 3) as usize);

    while xs >= (ys * 2) {
        mul_toom22(w_tmp, xp, ys, yp, ys, scratch);
        xs -= ys;
        xp = xp.offset(ys as isize);
        let cy = ll::add_n(wp, wp.as_const(), w_tmp.as_const(), ys);
        ll::copy_incr(w_tmp.offset(ys as isize).as_const(),
                      wp.offset(ys as isize),
                      ys);
        ll::incr(wp.offset(ys as isize), cy);

        wp = wp.offset(ys as isize);
    }

    if xs >= ys {
        mul_rec(w_tmp, xp, xs, yp, ys, scratch);
    } else {
        mul_rec(w_tmp, yp, ys, xp, xs, scratch);
    }

    let cy = ll::add_n(wp, wp.as_const(), w_tmp.as_const(), ys);
    ll::copy_incr(w_tmp.offset(ys as isize).as_const(),
                  wp.offset(ys as isize),
                  xs);
    ll::incr(wp.offset(ys as isize), cy);
}

/**
 * Squares the number in `{xp, xs}` storing the result in `{wp, xs*2}`.
 * This is slightly more efficient than regular multiplication with both
 * inputs the same.
 *
 * `{wp, xs*2}` must not overlap with `{xp, xs}`
 */
pub unsafe fn sqr(wp: LimbsMut, xp: Limbs, xs: i32) {
    debug_assert!(xs > 0);
    debug_assert!(!overlap(wp, 2*xs, xp, xs));

    if xs <= TOOM22_THRESHOLD {
        mul_basecase(wp, xp, xs, xp, xs);
    } else {
        let mut tmp = mem::TmpAllocator::new();
        let scratch = tmp.allocate((xs * 2) as usize);

        sqr_toom2(wp, xp, xs, scratch);
    }
}

#[inline(always)]
pub unsafe fn sqr_rec(wp: LimbsMut, xp: Limbs, xs: i32, scratch: LimbsMut) {
    if xs < TOOM22_THRESHOLD {
        mul_basecase(wp, xp, xs, xp, xs);
    } else {
        sqr_toom2(wp, xp, xs, scratch);
    }
}

unsafe fn sqr_toom2(wp: LimbsMut, xp: Limbs, xs: i32, scratch: LimbsMut) {
    // This is very similar to regular mul_toom22, however it is slightly more efficient
    // as it can take advantage of the coefficents being the same.
    //
    // Splitting x into x1, x0 to get x = x1*(B^n) + x0 means we get
    //
    //    x*x = (B^2n)*z2 + 2*(B^n)*z1 + z0
    //
    // Where:
    //   z0 = x0*x0
    //   z1 = x0*x1
    //   z2 = x1*x1

    let xh = xs >> 1;
    let xl = xs - xh;

    let x0 = xp;
    let x1 = xp.offset(xl as isize);

    let z0 = wp;
    let z1 = scratch;
    let z2 = wp.offset((xl * 2) as isize);
    let scratch_out = scratch.offset((xl * 2) as isize);

    // Calculate z1
    mul_rec(z1, x0, xl, x1, xh, scratch_out);
    // Calculate z0
    sqr_rec(z0, x0, xl, scratch_out);
    // Calculate z2
    sqr_rec(z2, x1, xh, scratch_out);

    // Calculate 2*z1
    let mut cy = ll::add_n(z1, z1.as_const(), z1.as_const(), xs);

    // wp now contains the result of (B^2n)*z2 + z0

    cy = cy + ll::add_n(wp.offset(xl as isize), wp.offset(xl as isize).as_const(), z1.as_const(), xs);

    ll::incr(wp.offset((xl + xs) as isize), cy);
}

#[cfg(test)]
mod test {

    #[test]
    fn test_mul_1() {
        use ll::limb_ptr::{Limbs, LimbsMut};
        use ll::limb::Limb;
        unsafe {
            let half_limb = 1 << (Limb::BITS-1);
            for &(a, l, x, x_c) in &[
                (&[1usize] as &[usize], 2, &[2usize] as &[usize], 0),
                (&[0, 1], 2, &[0, 2], 0),
                (&[1, 1], 2, &[2, 2], 0),
                (&[1, 1, 1], 2, &[2, 2, 2], 0),
                (&[1, 1, 1, 1], 2, &[2, 2, 2, 2], 0),
                (&[1, 2, 3, 4, 5], 2, &[2, 4, 6, 8, 10], 0),
                (&[half_limb], 2, &[0], 1),
                (&[0, half_limb], 2, &[0, 0], 1),
                (&[0, 0, half_limb], 2, &[0, 0, 0], 1),
                (&[0, 0, 0, half_limb], 2, &[0, 0, 0, 0], 1),
                (&[0, 0, 0, 0, half_limb], 2, &[0, 0, 0, 0, 0], 1),
                (&[!0], !0, &[1], !0-1),
            ] {
                let limbs = Limbs::new(a.as_ptr() as _, 0, a.len() as i32);
                let res_vec = vec!(0usize; a.len());
                let res = LimbsMut::new(res_vec.as_ptr() as _, 0, a.len() as i32);
                let Limb(carry) = super::mul_1(res, limbs, a.len() as _, Limb(l));
                assert_eq!(x_c, carry, "wrong carry testing {:?} * {}", a, l);
                assert_eq!(x, &*res_vec, "wrong result testing {:?} * {}", a, l);
            }
        }
    }

    /*
    #[test]
    fn test_addmul_2() {
        use super::addmul_2;
        use ll::limb_ptr::{Limbs, LimbsMut};
        use ll::limb::Limb;
        unsafe {
            for &(w, x, l1, l2, exp, exp_c) in &[
                (&[0usize, 0], &[0usize, 0], 0, 0, &[0usize, 0], (0, 0)),
                (&[1, 0], &[0, 0], 1, 0, &[1, 0], (0,0)),
                (&[0, 1], &[0, 0], 1, 0, &[0, 1], (0,0)),
                (&[0, 0], &[0, 0], 1, 0, &[0, 0], (0,0)),
                (&[0, 0], &[1, 0], 1, 0, &[1, 0], (0,0)),
                (&[1, 0], &[1, 0], 1, 0, &[2, 0], (0,0)),
                (&[0, 1], &[1, 0], 1, 0, &[1, 1], (0,0)),
                (&[0, 0], &[1, 0], 1, 0, &[1, 0], (0,0)),
                (&[0, 0], &[1, 0], 0, 1, &[0, 1], (0,0)),
                (&[0, 0], &[0, 1], 0, 1, &[0, 0], (1,0)),
                (&[0, 0], &[!0, !0], !0, !0, &[1, 0], (!0-1, !0)),
            ] {
                let x_vec = x.to_vec();
                let x_limbs = Limbs::new(x_vec.as_ptr() as _, 0, x.len() as i32);
                let w_vec = w.to_vec();
                let w_limbs = LimbsMut::new(w_vec.as_ptr() as _, 0, w_vec.len() as i32);
                let (Limb(c1), Limb(c2)) = addmul_2(w_limbs, x_limbs, x.len() as _, Limb(l1), Limb(l2));
                assert_eq!(exp, &*w_vec,
                           "wrong result testing {:?} += {:?} * {}:{}", w, x, l2, l1);
                assert_eq!(exp_c, (c1, c2),
                           "wrong carry testing {:?} += {:?} * {}:{}", w, x, l2, l1);
            }
        }
    }
    */

    #[test]
    fn test_mul_basecase() {
        use ll::limb_ptr::{Limbs, LimbsMut};
        unsafe {
            for &(x, y, exp) in &[
                (&[0usize, 0] as &[usize], &[0usize, 0] as &[usize], &[0usize, 0, 0, 0] as &[usize]),
                (&[1, 0], &[1, 0], &[1usize, 0, 0, 0]),
                (&[!0, !0], &[1, 0], &[!0, !0, 0, 0]),
                (&[!0, !0], &[!0, !0], &[1, 0, !0-1, !0]),
                (&[!0, !0, !0], &[!0, !0, !0], &[1, 0, 0, !0-1, !0, !0]),
                (&[1], &[1, 2, 3], &[1, 2, 3, 0]),
                (&[1], &[1, 2, 3, 4], &[1, 2, 3, 4, 0]),
                (&[0, 2], &[1, 2, 3, 4], &[0, 2, 4, 6, 8, 0]),
            ] {
                let x_vec = x.to_vec();
                let y_vec = y.to_vec();
                let w_vec = vec!(0usize; x.len()+y.len());
                let x_limbs = Limbs::new(x_vec.as_ptr() as _, 0, x.len() as i32);
                let y_limbs = Limbs::new(y_vec.as_ptr() as _, 0, y.len() as i32);
                let w_limbs = LimbsMut::new(w_vec.as_ptr() as _, 0, w_vec.len() as i32);
                super::mul_basecase(w_limbs, x_limbs, x.len() as _, y_limbs, y.len() as _);
                assert_eq!(exp, &*w_vec,
                           "wrong result testing {:?}*{:?}={:?} ", x, y, w_vec);
            }
        }
    }

    macro_rules! one_bench {
        ($size:expr, $name:ident, $what:expr) => {
            #[bench]
            fn $name(b: &mut ::test::Bencher) {
                use rand::Rng;
                use ll::limb::Limb;
                use ll::limb_ptr::{Limbs, LimbsMut};
                let mut rng = ::rand::thread_rng();
                unsafe {

                    let vx:Vec<Limb> = (0..$size).map(|_| Limb(rng.next_u64() as _)).collect();
                    let vy:Vec<Limb> = (0..$size).map(|_| Limb(rng.next_u64() as _)).collect();
                    let mut vz:Vec<Limb> = (0..2*$size).map(|_| Limb(rng.next_u64() as _)).collect();
                    let x = Limbs::new(vx.as_ptr(), 0, $size as i32);
                    let y = Limbs::new(vy.as_ptr(), 0, $size as i32);
                    let z = LimbsMut::new(vz.as_mut_ptr(), 0, 2*$size as i32);

                    b.iter(|| $what(z,x,$size,y));
                }
            }
        }
    }

    macro_rules! ladder { ($what:expr) => {
        one_bench!(2, size_0002, $what);
        one_bench!(4, size_0004, $what);
        one_bench!(8, size_0008, $what);
        one_bench!(16, size_0016, $what);
        one_bench!(32, size_0032, $what);
        one_bench!(64, size_0064, $what);
        one_bench!(128, size_0128, $what);
        one_bench!(256, size_0256, $what);
        one_bench!(512, size_0512, $what);
        one_bench!(1024, size_1024, $what);
        one_bench!(2048, size_2048, $what);
        one_bench!(4096, size_4096, $what);
        one_bench!(8192, size_8192, $what);
    }}

    mod mul_1 { ladder!(|z,x,xs,y:Limbs| super::super::mul_1(z, x, xs as i32, *y)); }
    mod addmul_1 { ladder!(|z,x,xs,y:Limbs| super::super::addmul_1(z, x, xs as i32, *y)); }
//    mod addmul_2 { ladder!(|z,x,xs,y:Limbs| super::super::addmul_2(z, x, xs as i32, *y, *y.offset(1))); }
    mod mul_basecase { ladder!(|z,x,xs,y| super::super::mul_basecase(z, x, xs as i32, y, xs as i32)); }
}
