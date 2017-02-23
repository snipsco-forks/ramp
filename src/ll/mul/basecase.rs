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

use ll::limb_ptr::{Limbs, LimbsMut};

/**
 * Multiplies `{xp, xs}` by `{yp, ys}`, storing the result to `{wp, xs + ys}`.
 *
 * `{wp, xs + ys}` must be disjoint from both inputs.
 */
pub unsafe fn mul_basecase(wp: LimbsMut, xp: Limbs, xs: i32, yp: Limbs, ys: i32) {
    if true {
        asm_x86_64(wp, xp, xs, yp, ys)
    } else {
        generic(wp, xp, xs, yp, ys)
    }
}

#[inline]
pub unsafe fn generic(mut wp: LimbsMut, xp: Limbs, xs: i32, mut yp: Limbs, mut ys: i32) {

    *wp.offset(xs as isize) = ::ll::mul_1(wp, xp, xs, *yp);
    wp = wp.offset(1);
    yp = yp.offset(1);
    ys -= 1;

    while ys > 0 {
        *wp.offset(xs as isize) = ::ll::addmul_1(wp, xp, xs, *yp);

        wp = wp.offset(1);
        yp = yp.offset(1);
        ys -= 1;
    }
}

#[inline]
#[allow(unused_assignments)]
pub unsafe fn asm_x86_64(wp: LimbsMut, xp: Limbs, xs: i32, yp: Limbs, ys: i32) {
    use ll::limb::Limb;
    let mut w:*mut _ = &mut *wp.offset(0);
    let mut x:*const _ = &*xp.offset(0);
    let mut y:*const _ = &*yp.offset(0);
    let mut xs = xs;
    let mut ys = ys;
    asm!("
                xorq %r11, %r11             // r11 <- 0 (overflow)
                mov (%r10), %rbx            // rbx <- yi
                mov %r8, %rcx               // rcx <- xs
            5:                              // initial mul loop top (rcx)
                mov (%rsi), %rax
                mulq %rbx
                add %r11, %rax
                adc $$0, %rdx
                mov %rax, (%rdi)
                mov %rdx, %r11

                add $$8, %rsi
                add $$8, %rdi
                dec %rcx
                jnz 5b                      // end of mul loop

                mov %r11, (%rdi)            // write final word

                movq %r8, %rax              // reset %rsi and %rdi
                shlq $$3, %rax
                subq %rax, %rsi
                subq %rax, %rdi

                addq $$8, %rdi              // offset rdi
                addq $$8, %r10              // offset r10
                dec %r9
                jz 6f

            4:                              // outer loop top (r9/ys)
                xorq %r11, %r11             // r11 <- 0 (overflow)
                mov %r8, %rcx               // rcx <- xs
                mov (%r10), %rbx            // rbx <- yi

                mov (%rsi), %rax
                mulq %rbx
                add %rax, (%rdi)
                adc $$0, %rdx
                mov %rdx, %r11

                add $$8, %rsi
                add $$8, %rdi
                dec %rcx
                jz 2f
            1:                              // inner addmul loop top (rcx)
                mov (%rsi), %rax
                mulq %rbx
                add %r11, %rax
                adc $$0, %rdx
                mov %rdx, %r11
                add %rax, (%rdi)
                adc $$0, %r11

                add $$8, %rsi
                add $$8, %rdi
                dec %rcx
                jnz 1b                      // end of addmul loop
            2:
                mov %r11, (%rdi)            // overwrite final word, no add

                movq %r8, %rax              // reset %rsi and %rdi
                shlq $$3, %rax
                subq %rax, %rsi
                subq %rax, %rdi

                addq $$8, %rdi              // offset rdi
                addq $$8, %r10              // offset r10
                dec %r9
                jnz 4b
            6:
        "
        : "=&{rdi}"(w), "=&{rsi}"(x), "=&{r8}"(xs), "=&{r10}"(y), "=&{r9}"(ys)
        : "0"(w), "1"(x), "2"(xs), "3"(y), "4"(ys)
        :   "rax", "rcx", "rdx", "rbx", "r11",
            "cc", "memory"
    );
}


#[cfg(test)]
mod test {
    macro_rules! t {
        ($func:ident) => {
            #[test]
            fn $func() {
                use ll::limb::Limb;
                use ll::limb_ptr::{Limbs, LimbsMut};
                unsafe {
                    let h:usize = 1<<(Limb::BITS-1);
                    for &(x, y, exp) in &[
                        (&[0usize, 0] as &[usize], &[0usize, 0] as &[usize], &[0usize, 0, 0, 0] as &[usize]),
                        (&[1, 0], &[1, 0], &[1usize, 0, 0, 0]),
                        (&[1], &[1], &[1, 0]),
                        (&[h], &[2], &[0, 1]),
                        (&[h, h], &[2], &[0, 1, 1]),
                        (&[2], &[h], &[0, 1]),
                        (&[2], &[h, h, h], &[0, 1, 1, 1]),
                        (&[1], &[1, 2, 3], &[1, 2, 3, 0]),
                        (&[1], &[1, 2, 3, 4], &[1, 2, 3, 4, 0]),
                        (&[1, 2, 3, 4], &[1], &[1, 2, 3, 4, 0]),
                        (&[0, 2], &[1, 2, 3, 4], &[0, 2, 4, 6, 8, 0]),
                        (&[!0, !0], &[1, 0], &[!0, !0, 0, 0]),
                        (&[!0, !0], &[!0, !0], &[1, 0, !0-1, !0]),
                        (&[!0, !0, !0], &[!0, !0, !0], &[1, 0, 0, !0-1, !0, !0]),
                    ] {
                        println!("{:?} x {:?}", x, y); 
                        let x_vec = x.to_vec();
                        let y_vec = y.to_vec();
                        let w_vec = vec!(0usize; x.len()+y.len());
                        let x_limbs = Limbs::new(x_vec.as_ptr() as _, 0, x.len() as i32);
                        let y_limbs = Limbs::new(y_vec.as_ptr() as _, 0, y.len() as i32);
                        let w_limbs = LimbsMut::new(w_vec.as_ptr() as _, 0, w_vec.len() as i32);
                        super::$func(w_limbs, x_limbs, x.len() as _, y_limbs, y.len() as _);
                        assert_eq!(exp, &*w_vec,
                                   "wrong result testing {:?}*{:?}={:?} ", x, y, w_vec);
                    }
                }
            }
        }
    }

    t!(generic);
    t!(asm_x86_64);
}
