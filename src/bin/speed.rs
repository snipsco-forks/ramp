#[macro_use]
#[cfg(feature="speed")]
extern crate clap;
extern crate rand;
extern crate framp as ramp;

#[cfg(not(feature="speed"))]
fn main() {
    panic!("ramp speed program can only be built with --features speed");
}


#[cfg(feature="speed")]
fn main() {
    run::run()
}

#[cfg(feature="speed")]
mod run {
    use clap;
    use std::collections::HashMap;
    use ramp::ll::limb_ptr::{Limbs, LimbsMut};

    macro_rules! routine_size_ns {
        ($hash:ident, $name: ident, $what:expr) => {
            fn $name(size:usize) -> u64 {
                use rand::Rng;
                use rand::SeedableRng;
                use ramp::ll::limb::Limb;
                use ramp::ll::limb_ptr::{Limbs, LimbsMut};
                use std::time::Instant;
                let mut rng = ::rand::StdRng::from_seed(&[12]);
                unsafe {
                    let vx:Vec<Limb> = (0..size).map(|_| Limb(rng.next_u64() as _)).collect();
                    let vy:Vec<Limb> = (0..size).map(|_| Limb(rng.next_u64() as _)).collect();
                    let mut vz:Vec<Limb> = (0..2*size).map(|_| Limb(rng.next_u64() as _)).collect();
                    let x = Limbs::new(vx.as_ptr(), 0, size as i32);
                    let y = Limbs::new(vy.as_ptr(), 0, size as i32 + 1);
                    let z = LimbsMut::new(vz.as_mut_ptr(), 0, 2*size as i32);

                    let mut tmp = ::ramp::mem::TmpAllocator::new();
                    let scratch = tmp.allocate((size * 2) as usize);

                    let start = Instant::now();
                    $what(z,x,size as i32,y,scratch);
                    let stop = Instant::now();
                    let elapsed = stop-start;
                    (elapsed.subsec_nanos() as u64 + elapsed.as_secs() * 1_000_000_000)
                }
            }
            $hash.insert(stringify!($name), $name);
        }

    }

    macro_rules! time_ns {
        ($n:expr, $what:block) => {{
            use std::time::Instant;
            let start = Instant::now();
            for _ in 0..$n {
                $what;
            }
            let stop = Instant::now();
            let elapsed = stop-start;
            (elapsed.subsec_nanos() as u64 + elapsed.as_secs() * 1_000_000_000) / $n as u64
        }}
    }

    macro_rules! map(
        { $($key:expr => $value:expr),+ } => { {
                let mut m = ::std::collections::HashMap::new();
                $( m.insert($key, $value); )+
                m
        } };
    );


    #[inline(never)]
    fn noop_non_inline(_:LimbsMut, _:Limbs, _:i32, _:Limbs, _:i32, _:LimbsMut ) { }

    pub fn run() {
        let mut routines_size: HashMap<&'static str, fn(usize) -> u64> = HashMap::new();

        {
        let mut h = &mut routines_size;
        routine_size_ns!(h, noop, {});
        routine_size_ns!(h, non_inline_call, |wp,xp,xn,yp,s| noop_non_inline(wp,xp,xn,yp,xn,s));

        routine_size_ns!(h, mul, |wp,xp,xn,yp,_| ::ramp::ll::mul::mul(wp,xp,xn,yp,xn));
        routine_size_ns!(h, mul_toom22, |wp,xp,xn,yp,s| ::ramp::ll::mul::mul_toom22(wp,xp,xn,yp,xn,s));
        routine_size_ns!(h, mul_basecase, |wp,xp,xn,yp,_| ::ramp::ll::mul::mul_basecase(wp,xp,xn,yp,xn));
        routine_size_ns!(h, mul_1, |wp,xp,xn,yp:Limbs,_| ::ramp::ll::mul::mul_1(wp,xp,xn,*yp));
        routine_size_ns!(h, addmul_1, |wp,xp,xn,yp:Limbs,_| ::ramp::ll::mul::addmul_1(wp,xp,xn,*yp));
        routine_size_ns!(h, addmul_1_generic, |wp,xp,xn,yp:Limbs,_| ::ramp::ll::mul::addmul_1::generic(wp,xp,xn,*yp));
        routine_size_ns!(h, addmul_1_generic_unroll_4, |wp,xp,xn,yp:Limbs,_| ::ramp::ll::mul::addmul_1::generic_unroll_4(wp,xp,xn,*yp));
        }

        let matches = clap_app!(speed =>
            (about: "ramp routine bench tool")
            (@arg SIZE: -s +required +takes_value "sizes to try. can be a range (1..10 or 1-10) or a list (1,2,4,8)")
            (@arg FACTOR: -f +takes_value "factor to step in range mode (defaults 1.1)")
            (@arg LOOPS: -l +takes_value "number of loops to run")
            (@arg ROUTINES: +multiple +required "routines to bench")
        ).get_matches();

        let loops = matches.value_of("LOOPS").unwrap().parse().unwrap();

        let routine = routines_size.get("noop").unwrap();
        let mut ns = 0;
        for _ in 0..loops {
            ns += routine(1);
        }
        let overhead = ns/loops;

        let sizes = sizes(&matches);
        for size in sizes {
            print!("{:10}", size);
            for r in matches.values_of("ROUTINES").unwrap() {
                let routine = routines_size.get(r).unwrap();
                let mut values = vec!();
                for _ in 0..loops {
                    values.push(routine(size));
                }
                values.sort();
                print!("\t{:10}", values[loops as usize/2].saturating_sub(overhead));
            }
            println!("");
        }
    }

    pub fn sizes(matches: &clap::ArgMatches) -> Vec<usize> {
        let size = matches.value_of("SIZE").unwrap().replace("..", "-");
        let factor: f32 = matches.value_of("FACTOR").unwrap_or("1.1").parse().unwrap();
        size.split(",")
            .flat_map(|f| if f.contains("-") {
                let mut beg_end = f.split("-");
                let mut i: usize = beg_end.next().unwrap().parse().unwrap();
                let end: usize = beg_end.next().unwrap().parse().unwrap();
                let mut points = vec![];
                while i <= end {
                    points.push(i);
                    i = ::std::cmp::max((i as f32 * factor) as usize, i + 1);
                }
                points.into_iter()
            } else {
                vec!(f.parse().unwrap()).into_iter()
            })
            .collect()
    }
}
