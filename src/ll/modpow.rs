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

use ll;
use ll::limb::Limb;
use mem;

use ll::limb_ptr::{Limbs, LimbsMut};

// w <- a^b [m] 
pub unsafe fn modpow_by_montgomery(wp:LimbsMut, r_limbs:i32, n:Limbs, nquote:Limbs, a:Limbs, bp:Limbs, bn: i32) {
    let k = 3;

    let mut tmp = mem::TmpAllocator::new();
    let scratch = tmp.allocate(2*r_limbs as usize); // for temp muls

    // base ^ 0..2^(k-1)
    let mut table = Vec::with_capacity(1 << k);
    let mut pow_0 = tmp.allocate(r_limbs as usize);
    *pow_0 = Limb(1);
    let pow_1 = tmp.allocate(r_limbs as usize);
    ll::copy_incr(a, pow_1, r_limbs as i32);
    table.push(pow_0);
    table.push(pow_1);
    for _ in 2..(1 << k) {
        let next = tmp.allocate(r_limbs as usize);
        {
            let previous = table.last().unwrap();
            montgomery_mul(next, r_limbs, pow_1.as_const(), previous.as_const(), n, nquote);
        }
        table.push(next);
    }

    let exp_bit_length = ll::base::num_base_digits(bp, bn, 2) as usize;
    let block_count = (exp_bit_length + k - 1) / k;
    for i in (0..block_count).rev() {
        let mut block_value: usize = 0;
        for j in 0..k {
            let p = i*k+j;
            if p < exp_bit_length && (*(bp.offset((p/Limb::BITS) as isize)) >> (p%Limb::BITS)) & Limb(1) == Limb(1) {
                block_value |= 1 << j;
            }
        }
        for _ in 0..k {
            montgomery_mul(scratch, r_limbs, wp.as_const(), wp.as_const(), n, nquote);
            ll::copy_incr(scratch.as_const(), wp, r_limbs);
        }
        if block_value != 0 {
            montgomery_mul(scratch, r_limbs, wp.as_const(), table[block_value].as_const(), n, nquote);
            ll::copy_incr(scratch.as_const(), wp, r_limbs);
        }
    }
}

// unsafe fn d(a:Limbs, s:i32) -> String{
//     (0..s).rev().map(|l| {
//         let Limb(limb) = *a.offset(l as isize);
//         if limb == 0 {
//             "".to_string()
//         } else {
//             format!(" {:x}", limb)
//         }
//     }).collect()
// }

unsafe fn montgomery_mul(wp:LimbsMut, r_limbs:i32, a:Limbs, b:Limbs, n:Limbs, nquote:Limbs) {
    let mut tmp = mem::TmpAllocator::new();
    let scratch_t = tmp.allocate(2*r_limbs as usize);
    let scratch_m = tmp.allocate(2*r_limbs as usize);
    let scratch_mn = tmp.allocate(2*r_limbs as usize);
    let scratch_x = tmp.allocate(2*r_limbs as usize);

    // t <- a*b
    ll::mul::mul(scratch_t, a, r_limbs, b, r_limbs);
    // a*b % R is a*b [r_limbs]

    // M <- (a*b % R) N'
    ll::mul::mul(scratch_m, scratch_t.as_const(), r_limbs, nquote, r_limbs);

    // MN <- M%R N
    ll::mul::mul(scratch_mn, scratch_m.as_const(), r_limbs, n, r_limbs);

    // X <- T+MN
    ll::addsub::add_n(scratch_x, scratch_t.as_const(), scratch_mn.as_const(), 2*r_limbs);

    if ll::cmp(scratch_x.as_const().offset(r_limbs as isize), n, r_limbs) != ::std::cmp::Ordering::Less {
        ll::addsub::sub_n(wp, scratch_x.offset(r_limbs as isize).as_const(), n, r_limbs);
    } else {
        ll::copy_incr(scratch_x.offset(r_limbs as isize).as_const(), wp, r_limbs);
    }
}

// w <- a^b [m]
pub unsafe fn modpow(mut wp:LimbsMut, mp:Limbs, mn:i32, ap:Limbs, an: i32, bp:Limbs, bn: i32) {
    let k = 7;

    let mut tmp = mem::TmpAllocator::new();
    let scratch = tmp.allocate(2*mn as usize); // for temp muls
    let scratch_q = tmp.allocate(mn as usize + 1); // for divrem quotient

    // base ^ 0..2^(k-1)
    let mut table = Vec::with_capacity(1 << k);
    let mut pow_0 = tmp.allocate(mn as usize);
    *pow_0 = Limb(1);
    let pow_1 = tmp.allocate(mn as usize);
    ll::copy_incr(ap, pow_1, an);
    table.push(pow_0);
    table.push(pow_1);
    for _ in 2..(1 << k) {
        let next = tmp.allocate(mn as usize);
        {
            let previous = table.last().unwrap();
            ll::mul::mul(scratch, pow_1.as_const(), mn, previous.as_const(), mn);
            ll::div::divrem(scratch_q, next, scratch.as_const(), 2*mn, mp, mn);
        }
        table.push(next);
    }

    *wp = Limb(1);
    let exp_bit_length = ll::base::num_base_digits(bp, bn, 2) as usize;
    let block_count = (exp_bit_length + k - 1) / k;
    for i in (0..block_count).rev() {
        let mut block_value: usize = 0;
        for j in 0..k {
            let p = i*k+j;
            if p < exp_bit_length && (*(bp.offset((p/Limb::BITS) as isize)) >> (p%Limb::BITS)) & Limb(1) == Limb(1) {
                block_value |= 1 << j;
            }
        }
        for _ in 0..k {
            ll::mul::sqr(scratch, wp.as_const(), mn);
            ll::div::divrem(scratch_q, wp, scratch.as_const(), 2*mn, mp, mn);
        }
        if block_value != 0 {
            ll::mul::mul(scratch, table[block_value].as_const(), mn, wp.as_const(), mn);
            ll::div::divrem(scratch_q, wp, scratch.as_const(), 2*mn, mp, mn);
        }
    }
}
