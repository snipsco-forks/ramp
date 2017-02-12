extern crate framp;
extern crate rand;
use rand::Rng;
use framp::ll::limb::Limb;
use framp::ll::limb_ptr::{Limbs, LimbsMut};

fn main() {
    let mut rng = ::rand::thread_rng();
    unsafe {

        let vx: Vec<Limb> = (0..1024).map(|_| Limb(rng.next_u64() as _)).collect();
        let vy: Vec<Limb> = (0..1024).map(|_| Limb(rng.next_u64() as _)).collect();
        let mut vz: Vec<Limb> = (0..2048).map(|_| Limb(rng.next_u64() as _)).collect();
        let x = Limbs::new(vx.as_ptr(), 0, 1024);
        let y = Limbs::new(vy.as_ptr(), 0, 1024);
        let z = LimbsMut::new(vz.as_mut_ptr(), 0, 2048);

        framp::ll::mul::mul_basecase(z, x, 1024, y, 1024);
    }
}
