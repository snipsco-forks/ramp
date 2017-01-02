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
pub unsafe fn modpow(mut wp:LimbsMut, mp:Limbs, mn:i32, ap:Limbs, an: i32, bp:Limbs, bn: i32) {
    let mut tmp = mem::TmpAllocator::new();

    // a mn-sized copy of ap,an
    let base = tmp.allocate(mn as usize);
    ll::copy_incr(ap, base, an);

    ll::copy_incr(ap, wp, an);

    let scratch = tmp.allocate(2*mn as usize);

    // this one is never used (quotient in modulo ops)
    let scratch_q = tmp.allocate(mn as usize + 1);

    let exp_bit_length = ll::base::num_base_digits(bp, bn, 2) as usize;

    for i in (0..exp_bit_length-1).rev() {
        ll::mul::sqr(scratch, wp.as_const(), mn);
        ll::div::divrem(scratch_q, wp, scratch.as_const(), 2*mn, mp, mn);
        if (*(bp.offset((i/Limb::BITS) as isize)) >> (i%Limb::BITS)) & Limb(1) == Limb(1) {
            ll::mul::mul(scratch, base.as_const(), mn, wp.as_const(), mn);
            ll::div::divrem(scratch_q, wp, scratch.as_const(), 2*mn, mp, mn);
        }
    }
}
