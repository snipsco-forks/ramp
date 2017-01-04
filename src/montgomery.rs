
use int::Int;

#[derive(Debug)]
pub struct Montgomery<'a> {
    pub modulus: &'a Int,
    pub modulus_quote: Int,
    pub bits: usize,
    pub lomask: Int,
    pub r: Int,
    pub r_inv: Int,
    pub tmp: Int,
}

impl<'a> Montgomery<'a> {
    #[allow(dead_code)]
    pub fn new(modulus: &'a Int) -> Montgomery<'a> {
        use ll::limb::Limb;
        let limbs_count = (modulus.bit_length() as usize + Limb::BITS - 1) / Limb::BITS;
        let bits = limbs_count * Limb::BITS;
        let r = Int::one() << bits;
        Montgomery {
            modulus_quote: &r - inverse(modulus.clone(), &r).unwrap(),
            modulus: modulus,
            bits: bits,
            lomask: r.clone()-1,
            r: r.clone(),
            r_inv: inverse(r, &modulus).unwrap(),
            tmp: Int::zero(),
        }
    }

    #[allow(dead_code)]
    fn redc(&mut self, a: &mut Int) {
        self.tmp.clone_from(&*a);
        *a *= &self.modulus_quote;
        *a &= &self.lomask;
        *a *= self.modulus;
        *a += &self.tmp;
        *a >>= self.bits;
        if (a as &Int) >= self.modulus {
            *a -= self.modulus
        }
    }

    #[allow(dead_code)]
    fn mul(&mut self, a: &mut Int, b: &Int) {
        *a *= b;
        self.redc(a);
    }

    #[allow(dead_code)]
    fn sq(&mut self, a: &mut Int) {
        self.tmp.clone_from(a);
        a.clone_from(&self.tmp.square());
        self.redc(a);
    }

    #[allow(dead_code)]
    pub fn to_montgomery(&self, a: &Int) -> Int {
        ((a * &self.r) % self.modulus)
    }

    #[allow(dead_code)]
    pub fn to_natural(&self, a: &Int) -> Int {
        a * &self.r_inv % self.modulus
    }
}

#[allow(dead_code)]
fn modpow_lr_k_ary_with_k_mgmy(base: &Int, exp: &Int, modulus: &Int, k: usize) -> Int {
    let mut m = Montgomery::new(modulus);

    // precompute base ^ 0..2^(k-1)
    let mut table = Vec::with_capacity(1 << k);
    table.push(Int::one());
    table.push(m.to_montgomery(base));
    for _ in 2..(1 << k) {
        let mut next = table.last().unwrap().clone();
        m.mul(&mut next, &table[1]);
        table.push(next);
    }

    let mut result = m.to_montgomery(&Int::one());
    let block_count = (exp.bit_length() as usize + k - 1) / k;
    for i in (0..block_count).rev() {
        let mut block_value: usize = 0;
        for j in 0..k {
            if exp.bit((i * k + j) as u32) {
                block_value |= 1 << j;
            }
        }
        for _ in 0..k {
            m.sq(&mut result);
        }

        if block_value != 0 {
            m.mul(&mut result, &table[block_value]);
        }
    }
    m.to_natural(&result)
}

#[cfg(test)]
mod test {

    use test::Bencher;

    static P: &'static str = "148677972634832330983979593310074301486537017973460461278300587514468301043894574906886127642530475786889672304776052879927627556769456140664043088700743909632312483413393134504352834240399191134336344285483935856491230340093391784574980688823380828143810804684752914935441384845195613674104960646037368551517";
    static Q: &'static str = "158741574437007245654463598139927898730476924736461654463975966787719309357536545869203069369466212089132653564188443272208127277664424448947476335413293018778018615899291704693105620242763173357203898195318179150836424196645745308205164116144020613415407736216097185962171301808761138424668335445923774195463";
    static N: &'static str = "446397596678771930935753654586920306936946621208913265356418844327220812727766442444894747633541329301877801861589929170469310562024276317335720389819531817915083642419664574530820516411614402061341540773621609718596217130180876113842466833544592377419546315874157443700724565446359813992789873047692473646165446397596678771930935753654586920306936946621208913265356418844327220812727766442444894747633541329301877801861589929170469310562045923774195463";
    static P_EXP_Q_MOD_N: &'static str = "167216127033575887543627597836645861047205125657210928573959751482755137615538210337351142826820586625192642801613712405599811895698660256697022034706036302526688254935463675298422321466416268553928456486375399618780536765018283218497477719051444372227826918812735482583824151162705395833327342668518742611088648794167631267226166034273943473474852640344643160320108048818901941885781997670470039501703327746459928325135708813764810716722046027109043738";

    use ::Int;
    fn parse_them() -> (Int, Int, Int, Int) {
        use std::str::FromStr;
        let p = Int::from_str(P).unwrap();
        let q = Int::from_str(Q).unwrap();
        let n = Int::from_str(N).unwrap();
        let x = Int::from_str(P_EXP_Q_MOD_N).unwrap();
        (p, q, n, x)
    }

    #[test]
    fn test_montgomery() {
        let (p, q, n, x) = parse_them();
        let mut m = super::Montgomery::new(&n);
        assert_eq!(m.to_natural(&m.to_montgomery(&p)), p);
        assert_eq!(m.to_natural(&m.to_montgomery(&q)), q);
        let mut mp = m.to_montgomery(&p);
        let mq = m.to_montgomery(&q);
        m.mul(&mut mp, &mq);
        assert_eq!(m.to_natural(&mp), (&p*&q) % &n);
        let mut mp = m.to_montgomery(&p);
        m.sq(&mut mp);
        assert_eq!(m.to_natural(&mp), (&p*&p) % &n);
        assert_eq!(super::modpow_lr_k_ary_with_k_mgmy(&p, &q, &n, 7), x);
    }

    #[test]
    fn test_montgomery_mul() {
        let n : Int = "1009".parse().unwrap();
        let a : Int = "2".parse().unwrap();
        let b : Int = "10".parse().unwrap();
        let x : Int = "15".parse().unwrap();
        assert_eq!(super::modpow_lr_k_ary_with_k_mgmy(&a, &b, &n, 3), x);
    }

}

#[allow(dead_code)]
struct GcdResult {
    gcd: Int,
    c1: Int,
    c2: Int,
}

#[allow(dead_code)]
fn extended_gcd(a: Int, b: Int) -> GcdResult {
    // Euclid's extended algorithm
    let (mut s, mut old_s) = (Int::zero(), Int::one());
    let (mut t, mut old_t) = (Int::one(), Int::zero());
    let (mut r, mut old_r) = (b, a);

    let mut tmp = Int::zero();
    while r != 0 {
        let quotient = &old_r / &r;
        tmp.clone_from(&r);
        r = &old_r - &quotient * r;
        old_r.clone_from(&tmp);
        tmp.clone_from(&s);
        s = &old_s - &quotient * s;
        old_s.clone_from(&tmp);
        tmp.clone_from(&t);
        t = &old_t - &quotient * t;
        old_t.clone_from(&tmp);
    }

    let _quotients = (t, s); // == (a, b) / gcd

    GcdResult {
        gcd: old_r,
        c1: old_s,
        c2: old_t,
    }
}

/// Find the standard representation of a (mod n).
fn normalize(a: Int, n: &Int) -> Int {
    let a = a % n;
    match a.cmp(&Int::zero()) {
        ::std::cmp::Ordering::Less => a + n,
        _ => a,
    }
}

/// Calculate the inverse of a (mod n).
fn inverse(a: Int, n: &Int) -> Option<Int> {
    let GcdResult { gcd, c1: c, c2: _ } = extended_gcd(a, n.clone());
    if gcd == 1 {
        Some(normalize(c, n))
    } else {
        None
    }
}

