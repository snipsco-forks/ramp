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
