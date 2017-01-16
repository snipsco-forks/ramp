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
pub unsafe fn modpow(wp: LimbsMut, r_limbs: i32, n: Limbs, nquote0: Limb, a: Limbs, bp: Limbs, bn: i32) {
    let k = 6;

    let mut tmp = mem::TmpAllocator::new();
    let t = tmp.allocate((2 * r_limbs + 1) as usize);
    let scratch_mul = tmp.allocate(2 * r_limbs as usize);

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
            mul(next,
                r_limbs,
                pow_1.as_const(),
                previous.as_const(),
                n,
                nquote0,
                t,
                scratch_mul);
        }
        table.push(next);
    }

    let exp_bit_length = ll::base::num_base_digits(bp, bn, 2) as usize;
    let block_count = (exp_bit_length + k - 1) / k;
    for i in (0..block_count).rev() {
        let mut block_value: usize = 0;
        for j in 0..k {
            let p = i * k + j;
            if p < exp_bit_length &&
               (*(bp.offset((p / Limb::BITS) as isize)) >> (p % Limb::BITS)) & Limb(1) == Limb(1) {
                block_value |= 1 << j;
            }
        }
        for _ in 0..k {
            sqr(wp, r_limbs, wp.as_const(), n, nquote0, t, scratch_mul);
        }
        if block_value != 0 {
            mul(wp,
                r_limbs,
                wp.as_const(),
                table[block_value].as_const(),
                n,
                nquote0,
                t,
                scratch_mul);
        }
    }
}

#[inline]
unsafe fn mul(wp: LimbsMut,
              r_limbs: i32,
              a: Limbs,
              b: Limbs,
              n: Limbs,
              nquote0: Limb,
              t: LimbsMut,
              scratch_mul: LimbsMut) {
    ll::mul::mul_rec(t, a, r_limbs, b, r_limbs, scratch_mul);
    redc(wp, r_limbs, n, nquote0, t)
}

#[inline]
unsafe fn sqr(wp: LimbsMut,
              r_limbs: i32,
              a: Limbs,
              n: Limbs,
              nquote0: Limb,
              t: LimbsMut,
              scratch_mul: LimbsMut) {
    ll::mul::sqr_rec(t, a, r_limbs, scratch_mul);
    redc(wp, r_limbs, n, nquote0, t)
}

#[inline]
pub unsafe fn redc(wp: LimbsMut, r_limbs: i32, n: Limbs, nquote0: Limb, t: LimbsMut) {
    let mut carry = 0;
    for i in 0..r_limbs {
        carry = 0;
        let m = (*t.offset(i as _)).0.wrapping_mul(nquote0.0 as _);
        for j in 0..r_limbs {
            let (h_mnj, l_mnj) = Limb(m).mul_hilo(*(n.offset(j as _)));
            let (s, c1) = t.offset((i + j) as _).add_overflow(l_mnj);
            let (s, c2) = s.add_overflow(Limb(carry));
            carry = c1 as ll::limb::BaseInt + c2 as ll::limb::BaseInt + h_mnj.0;
            *t.offset((i + j) as _) = s;
        }
        for j in (i + r_limbs)..(2 * r_limbs) {
            let (s, c) = t.offset(j as _).add_overflow(Limb(carry));
            carry = c as _;
            *t.offset(j as _) = s;
        }
    }
    if carry > 0 ||
       ll::cmp(t.offset(r_limbs as isize).as_const(), n, r_limbs) != ::std::cmp::Ordering::Less {
        ll::addsub::sub_n(wp, t.offset(r_limbs as isize).as_const(), n, r_limbs);
    } else {
        ll::copy_incr(t.offset(r_limbs as isize).as_const(), wp, r_limbs);
    }
}

pub fn inv1(x: Limb) -> Limb {
    let Limb(x) = x;
    let mut y = 1;
    for i in 2..(Limb::BITS) {
        if 1 << (i - 1) < (x.wrapping_mul(y) % (1 << i)) {
            y += 1 << i - 1;
        }
    }
    if 1 << (Limb::BITS - 1) < x.wrapping_mul(y) {
        y += 1 << Limb::BITS - 1;
    }
    Limb(y as _)
}

#[test]
fn test_inv1() {
    assert_eq!(inv1(Limb(23)).0.wrapping_mul(23), 1);
}

#[cfg(target_pointer_width = "64")]
#[test]
fn test_inv1_64() {
    assert_eq!(inv1(Limb(193514046488575)).0.wrapping_mul(193514046488575),
               1);
}
