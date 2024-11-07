mod curve;
mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;

pub use curve::*;
pub use engine::*;
pub use fq::*;
pub use fq12::*;
pub use fq2::*;
pub use fq6::*;
pub use fr::*;

pub const BN_X: u64 = 4965661367192848881;

// 6U+2 for in NAF form
pub const SIX_U_PLUS_2_NAF: [i8; 65] = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0,
    0, 1, 0, 1, 1,
];

#[cfg(test)]
mod tests {
    use crate::bn256::{multi_miller_loop, Fq2, Fq6, G1Affine, G2Affine, Gt, BN_X, G1, G2};
    use crate::{bn256, bn256::Fr};
    use bn256::fq::Fq;
    use bn256::fq12::Fq12;
    use ff::{Field, PrimeField};
    use group::prime::PrimeCurveAffine;
    use group::{Curve, Group};
    use num_bigint::{BigInt, BigUint, ToBigInt};
    use num_integer::{div_floor, Integer};
    use num_traits::real::Real;
    use num_traits::{FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
    use pairing::MillerLoopResult;
    use pasta_curves::arithmetic::CurveAffine;
    use rand_core::OsRng;
    use std::mem::swap;
    use std::ops::{Add, AddAssign, Div, Mul, Neg, Rem, Sub, SubAssign};
    use subtle::Choice;

    #[test]
    fn test_canonical_pairing_example() {
        let a_g1 = G1Affine::generator();
        let a_g2 = G2Affine::generator();

        let b_g1 = -G1Affine::generator();
        let b_g2 = G2Affine::generator();

        // e(a_g1, a_g2) == e(b_g1, b_g2)
        // Equivalently:
        // e(a_g1, a_g2) * e(-b_g1, b_g2) == 1
        let mut mml_res = multi_miller_loop(&[(&a_g1, &a_g2), (&b_g1, &b_g2)]);

        assert_eq!(Gt::identity(), mml_res.final_exponentiation());
    }


    // swaps
    fn fq2_swap(fq2: Fq2) -> Fq2 {
        Fq2::new(fq2.c1, fq2.c0)
    }

    fn fq6_swap(fq6: Fq6) -> Fq6 {
        Fq6::new(fq2_swap(fq6.c2), fq2_swap(fq6.c1), fq2_swap(fq6.c0))
    }

    fn fq12_swap(fq12: Fq12) -> Fq12 {
        Fq12::new(fq6_swap(fq12.c1), fq6_swap(fq12.c0))
    }

    fn g2_swap(g2: G2Affine) -> G2Affine {
        let mut g2_copy = g2.clone();
        g2_copy.x = fq2_swap(g2_copy.x);
        g2_copy.y = fq2_swap(g2_copy.y);
        g2_copy
    }



    // printers
    fn print_fq2(fq2: Fq2) {
        println!(
            "{:?}, {:?}",
            BigInt::from_str_radix(&format!("{:?}", fq2.c0).as_str()[2..], 16).unwrap(),
            BigInt::from_str_radix(&format!("{:?}", fq2.c1).as_str()[2..], 16).unwrap()
        );
    }

    fn print_fq6(fq6: Fq6) {
        print_fq2(fq6.c0);
        print_fq2(fq6.c1);
        print_fq2(fq6.c2);
    }

    fn print_fq12(fq12: Fq12) {
        print_fq6(fq12.c0);
        print_fq6(fq12.c1);
    }

    fn print_g1(g1: G1Affine) {
        println!(
            "{:?}, {:?}",
            BigInt::from_str_radix(&format!("{:?}", g1.x).as_str()[2..], 16).unwrap(),
            BigInt::from_str_radix(&format!("{:?}", g1.y).as_str()[2..], 16).unwrap()
        );
    }

    fn print_g2(g2: G2Affine) {
        print_fq2(g2.x);
        print_fq2(g2.y);
    }



    // curve parameters
    fn e() -> BigInt {
        BigInt::from_u128(29793968203157093288).unwrap()
    }

    fn px_x() -> BigInt {
        BigInt::from_str_radix(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583",
            10,
        )
            .unwrap()
    }

    fn rx_x() -> BigInt {
        BigInt::from_str_radix(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
            .unwrap()
    }

    fn hx_x() -> BigInt {
        BigInt::from_str_radix(
            "552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480",
            10,
        )
            .unwrap()
    }

    fn lambdax_x() -> BigInt {
        BigInt::from_str_radix(
            "10486551571378427818905133077457505975146652579011797175399169355881771981095211883813744499745558409789005132135496770941292989421431235276221147148858384772096778432243207188878598198850276842458913349817007302752534892127325269",
            10,
        )
            .unwrap()
    }

    fn mx_x() -> BigInt {
        BigInt::from_str_radix(
            "479095176016622842441988045216678740802490611077830385734835461768408136080710527913838926422556221801736192429704607292331370270163904428626529743670357",
            10,
        )
            .unwrap()
    }

    fn px_x_minus_1() -> BigInt {
        BigInt::from_str_radix(
            "21888242871839275222246405745257275088696311157297823662689037894645226208582",
            10,
        )
            .unwrap()
    }

    fn g1_sage() -> G1Affine {
        let mut out = G1Affine::identity();
        out.x = bn256::fq::Fq::from_str_vartime(
            "19491323635986486980056165026003970884581302300479364565163758691834883767296",
        )
            .unwrap();
        out.y = bn256::fq::Fq::from_str_vartime(
            "2503817206389621232991390790939417031444960302945150474681637705185779211401",
        )
            .unwrap();
        out
    }
    fn g2_sage() -> G2Affine {
        //let mut g2 = G2::identity();
        let x = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "11403269112307582471523194844678173363615200121780745962919905543513926078845",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "10022529265301880767558967801827554994678953177337994173174782310334418209951",
            )
                .unwrap(),
        );
        let y = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "7417909083002664933410862546938954664060641619680344911439335935535164894254",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "14403937293889182757621054345090826401263455856569175761852807173588543872656",
            )
                .unwrap(),
        );
        let z = Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::one());

        g2_force_affine(x, y, z)
    }


    // functionality
    fn to_naf(x_: BigInt) -> Vec<i32> {
        let mut x = x_.clone();
        let mut z = vec![];

        while x.ge(&BigInt::one()) {
            if (x.mod_floor(&BigInt::from(2u64))).eq(&BigInt::zero()) {
                z.push(BigInt::zero());
            } else {
                let mut zi = BigInt::from(2u64);

                zi.sub_assign(BigInt::from(x.mod_floor(&BigInt::from(4u64))));

                x.sub_assign(zi.clone());

                z.push(zi);
            }

            x = x.div_floor(&BigInt::from_u64(2).unwrap());
        }

        z.iter().map(|item| item.to_i32().unwrap()).collect()
    }

    fn g1_scalar_mul(g1: G1Affine, value: BigInt) -> G1Affine {
        let mut e = to_naf(value);
        e.reverse();
        e = e[1..].to_vec();

        let g1_initial = g1.to_curve();
        let mut result = g1_initial.clone();

        for (i, kb) in e.into_iter().enumerate() {
            result = result.double();
            if kb == 1 {
                result = result.add(&g1_initial);
            } else if kb == -1 {
                result = result.sub(&g1_initial);
            }
        }
        result.to_affine()
    }

    fn g2_scalar_mul(g2: G2Affine, value: BigInt) -> G2Affine {
        let mut e = to_naf(value);
        e.reverse();
        e = e[1..].to_vec();

        let g2_initial = g2.clone();
        let mut result = g2_initial.clone();

        for (i, kb) in e.into_iter().enumerate() {
            result = g2_double(result);
            if kb == 1 {
                result = g2_add(result, g2_initial);
            } else if kb == -1 {
                result = g2_add(result, g2_negate(g2_initial));
            }
        }

        result
    }

    fn g2_force_affine(x: Fq2, y: Fq2, z: Fq2) -> G2Affine {
        if z == Fq2::zero()
            || z == Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::one())
            || z == Fq2::new(bn256::fq::Fq::one(), bn256::fq::Fq::zero())
        {
            let mut result = G2Affine::identity();
            result.x = x;
            result.y = y;
            return result;
        }

        let zinv = fq2_inverse(z);
        let zinv2 = fq2_mul(zinv, zinv);
        let zinv3 = fq2_mul(zinv2, zinv);

        let x_out = fq2_mul(x, zinv2);
        let y_out = fq2_mul(y, zinv3);

        let mut out = G2Affine::identity();
        out.x = x_out;
        out.y = y_out;
        out
    }

    fn g2_double(a: G2Affine) -> G2Affine {
        // http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

        let A = fq2_mul(a.x, a.x);
        let B = fq2_mul(a.y, a.y);
        let C = fq2_mul(B, B);

        let t = a.x + B;
        let t = fq2_mul(t, t);

        let D = t - A;
        let D = D - C;
        let D = D.double();
        let E = A.double();
        let E = E + A;
        let F = fq2_mul(E, E);

        let C8 = C.double().double().double();

        let c_x = F - D.double();
        let c_y = D - c_x;
        let c_y = fq2_mul(E, c_y);
        let c_y = c_y - C8;
        let c_z = fq2_mul(a.y, Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::one()));
        let c_z = c_z.double();

        g2_force_affine(c_x, c_y, c_z)
    }

    fn g2_negate(a: G2Affine) -> G2Affine {
        let mut out = a.clone();
        out.y = fq2_negative_of(out.y);
        out
    }

    fn g2_add(a: G2Affine, b: G2Affine) -> G2Affine {
        // http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl

        // TODO: handle infinite point cases

        let a_z = Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::one());
        let b_z = Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::one());

        let z1z1 = fq2_mul(a_z, a_z);
        let z2z2 = fq2_mul(b_z, b_z);

        let u1 = fq2_mul(z2z2, a.x);
        let u2 = fq2_mul(z1z1, b.x);
        let h = u2 - u1;

        let s1 = fq2_mul(a.y, b_z);
        let s1 = fq2_mul(s1, z2z2);
        let s2 = fq2_mul(b.y, a_z);
        let s2 = fq2_mul(s2, z1z1);
        let r = s2 - s1;

        if (h.is_zero().unwrap_u8() == 1u8) && (r.is_zero().unwrap_u8() == 1u8) {
            return g2_double(a);
        }

        let r = r.double();
        let i = fq2_mul(h, h);
        let i = i.double().double();
        let j = fq2_mul(h, i);
        let V = fq2_mul(u1, i);

        let c_x = fq2_mul(r, r);
        let c_x = c_x - j;
        let c_x = c_x - V.double();

        let tmp = fq2_mul(s1, j);
        let tmp = tmp.double();
        let c_y = V - c_x;
        let c_y = fq2_mul(r, c_y);
        let c_y = c_y - tmp;

        let c_z = a_z + b_z;
        let c_z = fq2_mul(c_z, c_z);
        let c_z = c_z - z1z1;
        let c_z = c_z - z2z2;
        let c_z = fq2_mul(c_z, h);
        g2_force_affine(c_x, c_y, c_z)
    }

    /*
    fn fq2_exp(fq2: Fq2, k: BigInt) -> Fq2 {
        let mut e = to_naf(k);
        e.reverse();
        e = e[1..].to_vec();

        let mut R = fq2.clone();
        for kb in e {
            R.square_assign();
            if kb == 1i32 {
                R *= R;
            } else if kb == -1i32 {
                R *= R.invert().unwrap();
            }
        }
        R
    }*/

    fn line_add(t: G2Affine, p: G2Affine) -> (Fq2, Fq2) {
        let x1 = t.x;
        let y1 = t.y;

        let x2 = p.x;
        let y2 = p.y;

        // slope: alpha = (y2 - y1) / (x2 - x1)
        let mut x2_x1 = fq2_inverse(x2 - x1);
        let mut y2_y1 = y2 - y1;
        let mut alpha = fq2_mul(y2_y1, x2_x1);
        let bias = y1 - fq2_mul(alpha, x1);
        (alpha, bias)
    }

    fn line_double(inner: G2Affine) -> (Fq2, Fq2) {
        let mut x_initial = inner.x.clone();
        let mut y_initial = inner.y.clone();
        let mut x = inner.x.clone();
        let mut y = inner.y.clone();

        x.square_assign();
        x *= Fq2::from(3);
        y *= Fq2::from(2);
        y = y.invert().unwrap();

        let alpha = x * y;
        let tmp = fq2_mul(alpha.clone(), x_initial);
        let bias = y_initial - tmp;

        (alpha, bias)
    }

    fn line_function(g2: G2Affine) -> Vec<(Fq2, Fq2)> {
        // TODO: simplify
        let mut naf_digits = to_naf(e());
        naf_digits.reverse();
        naf_digits = naf_digits[1..].to_vec();

        let mut L = vec![];

        let initial = g2.clone();
        let mut T = initial.clone();

        // TODO: optimize to_curve / to_affine invocations on T
        for (i, digit) in naf_digits.into_iter().enumerate() {
            let (alpha, bias) = line_double(T);
            T = g2_double(T);
            L.push((alpha, bias));

            if digit * digit == 1 {
                let qt = if digit == 1 { initial } else { -initial };
                let (alpha, bias) = line_add(T, qt);
                T = g2_add(T, qt);
                L.push((alpha, bias));
            }
        }
        assert_eq!(T, g2_scalar_mul(initial, e()));

        let mut pi_1_q = g2_scalar_mul(g2, px_x());
        let mut pi_2_q = g2_scalar_mul(g2, px_x().pow(2));
        let mut pi_3_q = g2_scalar_mul(g2, px_x().pow(3));

        // Here we use is_on_curve() method from halo2curves, so we need to swap c0 and c1 in order
        // to get halo2curves representation of the G2Affine point from sage
        assert_eq!(g2_swap(pi_1_q).is_on_curve().unwrap_u8(), 1u8);
        assert_eq!(g2_swap(pi_2_q).is_on_curve().unwrap_u8(), 1u8);
        assert_eq!(g2_swap(pi_3_q).is_on_curve().unwrap_u8(), 1u8);

        let (alpha, bias) = line_add(T, pi_1_q);
        T = g2_add(T, pi_1_q);
        L.push((alpha, bias));
        assert_eq!(T, g2_scalar_mul(initial, e() + px_x()));

        let (alpha, bias) = line_add(T, -pi_2_q);
        T = g2_add(T, -pi_2_q);
        L.push((alpha, bias));

        let k = e() + px_x() - px_x().pow(2);

        if k.clone().is_positive() {
            assert_eq!(T, g2_scalar_mul(initial, k.clone()));
        }

        if k.clone().is_negative() || k.clone().is_zero() {
            assert_eq!(
                T,
                g2_scalar_mul(initial, rx_x() - ((-k.clone()).mod_floor(&rx_x())))
            );
        }

        let alpha = Fq2::zero();
        // FIXME: double check if this is correct
        let bias = T.x;

        L.push((alpha, bias));
        L
    }

    fn fq2_mul(a: Fq2, b: Fq2) -> Fq2 {
        // If we want we can use halo2curves implementation alternatively:
        // let tmp = fq2_swap(fq2_swap(alpha.clone()) * fq2_swap(x_initial));
        let vy = a.c1 * b.c1;
        let vx = a.c0 * b.c0;
        let c0 = vy - vx;
        let c1 = (a.c0 + a.c1) * (b.c0 + b.c1) - vy - vx;
        Fq2::new(c1, c0)
    }

    fn fq_additive_inverse(a: bn256::fq::Fq) -> bn256::fq::Fq {
        let mod_1 =
            bn256::fq::Fq::from_str_vartime(px_x_minus_1().to_str_radix(10).as_str()).unwrap();
        let mut output = mod_1 - a;
        output += bn256::fq::Fq::one();
        output
    }

    fn fq2_inverse(a: Fq2) -> Fq2 {
        let t = a.c0.square() + a.c1.square();
        let inv = t.invert().unwrap();
        let c_x = fq_additive_inverse(a.c0) * inv;
        let c_y = a.c1 * inv;
        Fq2::new(c_x, c_y)
    }

    fn line_evaluation(alpha: Fq2, bias: Fq2, P: G1Affine) -> (Fq2, Fq2, Fq2) {
        let one = fq2_negative_of(bias);
        let two = fq2_negative_of(alpha);
        let two = fq2_mul_scalar(two, P.x);
        let three = Fq2::new(bn256::fq::Fq::zero(), P.y);
        (one, two, three)
    }

    fn fq2_negative_of(a: Fq2) -> Fq2 {
        let mut a_copy = a.clone();
        a_copy.c0 = fq_additive_inverse(a.c0);
        a_copy.c1 = fq_additive_inverse(a.c1);
        a_copy
    }

    fn fq2_mul_scalar(a: Fq2, scalar: bn256::fq::Fq) -> Fq2 {
        let mut a_copy = a.clone();
        a_copy.c0 = a.c0 * scalar;
        a_copy.c1 = a.c1 * scalar;
        a_copy
    }

    fn fq2_mul_beta(a: Fq2) -> Fq2 {
        // (beta + y)(3 + i) = 3beta + 3y - x + yi = (3x + y)i + (3y - x)
        let nine = a.c0 + a.c0 + a.c0 + a.c0 + a.c0 + a.c0 + a.c0 + a.c0 + a.c0;
        let tx = nine + a.c1;
        let nine = a.c1 + a.c1 + a.c1 + a.c1 + a.c1 + a.c1 + a.c1 + a.c1 + a.c1;
        let ty = nine - a.c0;
        Fq2::new(tx, ty)
    }

    fn fq6_mul_tau(a: Fq6) -> Fq6 {
        let tx = a.c1;
        let ty = a.c2;
        let tz = fq2_mul_beta(a.c0);
        Fq6::new(tx, ty, tz)
    }

    fn fq6_mul(a: Fq6, b: Fq6) -> Fq6 {
        // Algorithm 13 from http://eprint.iacr.org/2010/354.pdf

        /* TODO
        if self.x.is_zero():
            if self.y.is_zero():
                return b.mul_scalar(self.z)

            t0 = (b.z * self.z)
            t1 = (b.y * self.y)

            # tz = (b.x + self.y) * (self.y)
            tz = (b.x + b.y) * (self.y)
            tz -= t1
            tz = tz.mul_beta()
            tz += t0

            ty = (b.y + b.z) * (self.y + self.z)
            ty -= t0
            ty -= t1

            tx = (b.x) * (self.z)
            tx += t1

            return Fp6(tx, ty, tz)

        if b.x.is_zero():
            if b.y.is_zero():
                return self.mul_scalar(b.z)

            t0 = (self.z * b.z)
            t1 = (self.y * b.y)

            tz = (self.x + self.y) * (b.y)
            tz -= t1
            tz = tz.mul_beta()
            tz += t0

            ty = (self.y + self.z) * (b.y + b.z)
            ty -= t0
            ty -= t1

            tx = (self.x) * (b.z)
            tx += t1

            return Fp6(tx, ty, tz)
        */

        let t0 = fq2_mul(a.c2, b.c2);
        let t1 = fq2_mul(a.c1, b.c1);
        let t2 = fq2_mul(a.c0, b.c0);

        let mut tz = fq2_mul(a.c0 + a.c1, b.c0 + b.c1);
        tz -= t1;
        tz -= t2;
        tz = fq2_mul_beta(tz);
        tz += t0;

        let mut ty = fq2_mul(a.c1 + a.c2, b.c1 + b.c2);
        ty -= t0;
        ty -= t1;
        ty += fq2_mul_beta(t2);

        let mut tx = fq2_mul(a.c0 + a.c2, b.c0 + b.c2);
        tx -= t0;
        tx += t1;
        tx -= t2;

        Fq6::new(tx, ty, tz)
    }

    fn fq12_mul(a: Fq12, b: Fq12) -> Fq12 {
        // Karatsuba
        let axbx = fq6_mul(a.c0, b.c0);
        let axby = fq6_mul(a.c0, b.c1);
        let aybx = fq6_mul(a.c1, b.c0);
        let ayby = fq6_mul(a.c1, b.c1);

        Fq12::new(axby + aybx, ayby + fq6_mul_tau(axbx))
    }

    fn fq12_mul_line_base(r: Fq12, abc: (Fq2, Fq2, Fq2)) -> Fq12 {
        let fl = Fq12::new(
            Fq6::new(Fq2::zero(), abc.0, abc.1),
            Fq6::new(Fq2::zero(), Fq2::zero(), abc.2),
        );

        fq12_mul(r, fl)
    }

    fn multi_miller_loop_sage(terms: &[(&G1Affine, Vec<(Fq2, Fq2)>)]) -> Fq12 {
        let mut lc = 0;
        let mut f = Fq12::new(
            Fq6::new(
                Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
            ),
            Fq6::new(
                Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::one()),
            ),
        );

        let mut naf_digits = to_naf(e());
        naf_digits.reverse();
        naf_digits = naf_digits[1..].to_vec();

        for (i, digit) in naf_digits.into_iter().enumerate() {
            f = fq12_mul(f, f);

            for (j, (P, L)) in terms.into_iter().enumerate() {
                let (alpha, bias) = L[lc];
                let le = line_evaluation(alpha, bias, **P);
                f = fq12_mul_line_base(f, le);
                //f_list.push(f);

                if digit * digit == 1 {
                    let (alpha, bias) = L[lc + 1];
                    let le = line_evaluation(alpha, bias, **P);
                    f = fq12_mul_line_base(f, le);
                    //f_list.push(f);
                }
            }

            if digit == 0 {
                lc += 1;
            } else {
                lc += 2;
            }
        }

        // frobenius part
        for (j, (P, L)) in terms.into_iter().enumerate() {
            for k in 0..3 {
                let (alpha, bias) = L[lc + k];
                if k == 2 {
                    let eval = Fq12::new(
                        Fq6::new(
                            Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                            Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                            Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                        ),
                        Fq6::new(
                            Fq2::new(bn256::fq::Fq::zero(), bn256::fq::Fq::zero()),
                            fq2_negative_of(bias),
                            Fq2::new(bn256::fq::Fq::zero(), (**P).x),
                        ),
                    );
                    f = fq12_mul(f, eval);
                } else {
                    let le = line_evaluation(alpha, bias, **P);
                    f = fq12_mul_line_base(f, le);
                    //f_list.push(f);
                }
            }
        }

        lc += 3;
        assert_eq!(lc, terms[0].1.len());

        f
    }

    #[test]
    fn test_development() {
        let mut P1_sage = G1Affine::identity();
        P1_sage.x = bn256::fq::Fq::from_str_vartime(
            "17436239228436682145589725535260983283023841756095494372001414785830786557350",
        )
            .unwrap();
        P1_sage.y = bn256::fq::Fq::from_str_vartime(
            "2662552661475325419730296320054178759925559877515598138032134890479686302160",
        )
            .unwrap();

        let mut P2_sage = G1Affine::identity();
        P2_sage.x = bn256::fq::Fq::from_str_vartime(
            "17889629836789884053832158472344387789470457240319874813021323254969957756159",
        )
            .unwrap();
        P2_sage.y = bn256::fq::Fq::from_str_vartime(
            "3257012016312567276892866616568971437183359845703546812812284257267928868651",
        )
            .unwrap();

        assert_eq!(P1_sage, g1_scalar_mul(g1_sage(), BigInt::from(10i32)));
        assert_eq!(P2_sage, g1_scalar_mul(g1_sage(), BigInt::from(20i32)));

        let mut Q1_sage = G2Affine::identity();
        Q1_sage.x = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "2436744708922661243601109620162217258802612295655460342568049867468349557571",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "8467960934623505539437205290455186450587694510185411693257898784957175626315",
            )
                .unwrap(),
        );
        Q1_sage.y = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "17485882210820201413430435068015453483097088005813491924936525801298239061931",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "13486534423526602372174565913947995051694049622729666749442854723786326188044",
            )
                .unwrap(),
        );
        assert_eq!(Q1_sage, g2_scalar_mul(g2_sage(), BigInt::from(20i32)));

        let mut Q2_sage = G2Affine::identity();
        Q2_sage.x = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "18035812090211841711412088385265555624419102059185285619014957854927721278967",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "21723319147137606129563825226097395692782138078544663621240945871371724147820",
            )
                .unwrap(),
        );
        Q2_sage.y = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "554826965740810442399515967370827135896673026540856695059970955650940244188",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "10614327209813758539264379722287304292612535625145217656615353889802477615341",
            )
                .unwrap(),
        );
        assert_eq!(
            Q2_sage,
            g2_negate(g2_scalar_mul(g2_sage(), BigInt::from(10i32)))
        );

        let actual = line_add(Q1_sage, Q2_sage);
        assert_eq!(
            (
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("9898601632308850422125980073089659018033565770550791548603707476086123314398").unwrap(),
                    bn256::fq::Fq::from_str_vartime("14656042077947547196455876739877254672548486917221347731705635555435456288968").unwrap(),
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("19639791820696888430053130883591037135690871699014898409169105874022604989210").unwrap(),
                    bn256::fq::Fq::from_str_vartime("8759527323896261661677490200020184807605654202534521249636923264228006710187").unwrap()
                )
            ), actual);

        let actual = line_double(Q1_sage);
        assert_eq!(
            (
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("6217129506805939686487959038568043494839607178290145446577698642060725470602").unwrap(),
                    bn256::fq::Fq::from_str_vartime("19186596770352736118402211078625824463146104199773832148960967097929928937906").unwrap(),
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("16565183222735183739989499570062744955448770039312045718948385333082241858530").unwrap(),
                    bn256::fq::Fq::from_str_vartime("10683666479235621273392625200688677894541755144918312196152722001435717044713").unwrap()
                )
            ), actual);

        let actual = g2_double(Q2_sage);
        let mut expected = G2Affine::identity();
        expected.x = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "2436744708922661243601109620162217258802612295655460342568049867468349557571",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "8467960934623505539437205290455186450587694510185411693257898784957175626315",
            )
                .unwrap(),
        );
        expected.y = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "4402360661019073808815970677241821605599223151484331737752512093346987146652",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "8401708448312672850071839831309280037002261534568156913246183170858900020539",
            )
                .unwrap(),
        );

        assert_eq!(actual, expected);

        let actual = g2_add(Q1_sage, Q2_sage);
        let mut expected = G2Affine::identity();
        expected.x = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "18035812090211841711412088385265555624419102059185285619014957854927721278967",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "21723319147137606129563825226097395692782138078544663621240945871371724147820",
            )
                .unwrap(),
        );
        expected.y = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "21333415906098464779846889777886447952799638130756966967629066938994285964395",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "11273915662025516682982026022969970796083775532152606006073684004842748593242",
            )
                .unwrap(),
        );

        assert_eq!(actual, expected);

        let v = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "7417909083002664933410862546938954664060641619680344911439335935535164894254",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "14403937293889182757621054345090826401263455856569175761852807173588543872656",
            )
                .unwrap(),
        );

        let actual = fq2_negative_of(v);
        assert_eq!(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime(
                    "14470333788836610288835543198318320424635669537617478751249701959110061314329"
                )
                    .unwrap(),
                bn256::fq::Fq::from_str_vartime(
                    "7484305577950092464625351400166448687432855300728647900836230721056682335927"
                )
                    .unwrap()
            ),
            actual
        );

        let scalar = bn256::fq::Fq::from_str_vartime(
            "11403269112307582471523194844678173363615200121780745962919905543513926078845",
        )
            .unwrap();
        let actual = fq2_mul_scalar(v, scalar);
        assert_eq!(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime(
                    "7777384869092974241151079234238042410102508151892670643296383728229952533341"
                )
                    .unwrap(),
                bn256::fq::Fq::from_str_vartime(
                    "15740278586851715236335832894897416588509611562745531279098739159622574985629"
                )
                    .unwrap()
            ),
            actual
        );

        let actual = fq2_mul_beta(v);
        assert_eq!(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime(
                    "15500390425395341491579600031769593111720296961798808976739716909469349295193"
                )
                    .unwrap(),
                bn256::fq::Fq::from_str_vartime(
                    "12776312202803603773946597832592107503828905302953118631790739153535598916735"
                )
                    .unwrap()
            ),
            actual
        );

        let c0 = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "10248533281718687949696824861893906153006374018662717467375012751458360998749",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "16268156138197650612220956241172424778237127042309366985643003788250282721486",
            )
                .unwrap(),
        );
        let c1 = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "10125742394723339129880904802578269041133121936479011800291560529806659999223",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "1001675692397768398752533825968817244277481576409312917938900052133402393721",
            )
                .unwrap(),
        );
        let c2 = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "609610916123020890306391684018490082742641761110335093396493730048749660486",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "6090040222519325545030710964562589156473974467251038353775113433346320413849",
            )
                .unwrap(),
        );
        let c3 = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "21487194175735143898196256155845945983665814274918516790390963118326848922147",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "19963123303780593155840206717811652868087455511364366978073666602467509690418",
            )
                .unwrap(),
        );
        let c4 = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "5181263637101988486167218786997919676041719614354844817127017076758935925670",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "3431319842227613109421121927117318108215095992850705990130037913077038368788",
            )
                .unwrap(),
        );
        let c5 = Fq2::new(
            bn256::fq::Fq::from_str_vartime(
                "4025620765675370512048013633582632930081904187916626943601242132348839252752",
            )
                .unwrap(),
            bn256::fq::Fq::from_str_vartime(
                "5395205492065411747182189390509072658561418135451638572552184176878331010255",
            )
                .unwrap(),
        );

        let fq6_c0 = Fq6::new(c0, c1, c2);
        let fq6_c1 = Fq6::new(c3, c4, c5);
        let fq12 = Fq12::new(fq6_c0, fq6_c1);

        let actual = fq6_mul_tau(fq6_c0);
        assert_eq!(Fq6::new(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("10125742394723339129880904802578269041133121936479011800291560529806659999223").unwrap(),
                bn256::fq::Fq::from_str_vartime("1001675692397768398752533825968817244277481576409312917938900052133402393721").unwrap()
            ),
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("609610916123020890306391684018490082742641761110335093396493730048749660486").unwrap(),
                bn256::fq::Fq::from_str_vartime("6090040222519325545030710964562589156473974467251038353775113433346320413849").unwrap()
            ),
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("20951984186308741270506757017188479800509248581082529541261966972794626875895").unwrap(),
                bn256::fq::Fq::from_str_vartime("4835414731024516226813346837114266318949902418334643427277793974922826243127").unwrap()
            )
        ), actual);

        let actual = fq6_mul(fq6_c0, fq6_c1);
        assert_eq!(Fq6::new(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("17963686575674329183376140185280257330529413475832137100945531006693560344090").unwrap(),
                bn256::fq::Fq::from_str_vartime("17103024242318889364474531034832282411402848036059875219701172051391279864751").unwrap()
            ),
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("20635576593895534978130931469494361992329755594342296180882287102242485389595").unwrap(),
                bn256::fq::Fq::from_str_vartime("14687559097990452022801626604269559769630849304339150243143389297733541314395").unwrap()
            ),
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("18589950296866907975432879881736295960775030528683008526839918053376353452497").unwrap(),
                bn256::fq::Fq::from_str_vartime("1247725565779871246730782119434657049250543609736130888440159325440707230418").unwrap()
            )
        ), actual);

        let actual = fq12_mul(fq12, fq12);
        assert_eq!(Fq12::new(
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("14039130279509383144505874625303239572362515794366450539202024118741894479597").unwrap(),
                    bn256::fq::Fq::from_str_vartime("12317805612798503506702656324407289734109384914821926776713306208137333520919").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("19382910315951794734015457193731448895963200031386768699075536309839744570607").unwrap(),
                    bn256::fq::Fq::from_str_vartime("7486875324141628823356847463281844450565387451380476823597740700821856420207").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("15291657721894540728619354018215316832853749900068193390990798212107480696411").unwrap(),
                    bn256::fq::Fq::from_str_vartime("2495451131559742493461564238869314098501087219472261776880318650881414460836").unwrap()
                )
            ),
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("2949848866197155337881871387685395926098030206452722189454338691805468917875").unwrap(),
                    bn256::fq::Fq::from_str_vartime("3553704550311135302611682475262530596908939945896080877811304118934030306065").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("7522026249227091840777417884004471979183599748088190623943070757549008608516").unwrap(),
                    bn256::fq::Fq::from_str_vartime("4752431194108764974438947165608648388085578076667262698606995265022072055442").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("20976210973526535106000218646618966986561326314129333588288420626310888577751").unwrap(),
                    bn256::fq::Fq::from_str_vartime("6953137387553377021869901766560127966060787179798162717160245539851491987429").unwrap()
                )
            )
        ), actual);

        let actual = fq12_mul_line_base(fq12, (
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("7777384869092974241151079234238042410102508151892670643296383728229952533341").unwrap(),
                bn256::fq::Fq::from_str_vartime("15740278586851715236335832894897416588509611562745531279098739159622574985629").unwrap()
            ),
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("20976210973526535106000218646618966986561326314129333588288420626310888577751").unwrap(),
                bn256::fq::Fq::from_str_vartime("6953137387553377021869901766560127966060787179798162717160245539851491987429").unwrap()
            ),
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("14039130279509383144505874625303239572362515794366450539202024118741894479597").unwrap(),
                bn256::fq::Fq::from_str_vartime("12317805612798503506702656324407289734109384914821926776713306208137333520919").unwrap()
            )
        )
        );

        assert_eq!(Fq12::new(
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("14497372640154767938611099833932970340425757116277838117759615948094505211665").unwrap(),
                    bn256::fq::Fq::from_str_vartime("10282685723694816319133989374244602483859080709184087264567673833398995527741").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("8781796806787585092246316507665894587955051782828094008160738452255453821474").unwrap(),
                    bn256::fq::Fq::from_str_vartime("2055002370094466055406819187745197444209903848338323925440228384480359919559").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("20942500697119399512958070668451460092846871076311750451623959845078976794345").unwrap(),
                    bn256::fq::Fq::from_str_vartime("5344384980329075430872766168033161026834135066859716306510980180661384227096").unwrap()
                )
            ),
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("19149574538957921663043171550126182712315163359548676228444190507901769535535").unwrap(),
                    bn256::fq::Fq::from_str_vartime("2710262164390705464425409906656712810493705348412593118632772003135089084334").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("2775273036634849119074496113823008877208633800235788358822778359595285799915").unwrap(),
                    bn256::fq::Fq::from_str_vartime("20534417807606300188897454694864813727932624862189471978801381713934112032707").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("780544543800247237028432553069626896397197832944391721928747407672348123521").unwrap(),
                    bn256::fq::Fq::from_str_vartime("2220683654809138416095811755684847401448955848710261692900407609660774124105").unwrap()
                )
            )
        ), actual);

        let L1 = line_function(Q1_sage);
        let L2 = line_function(Q2_sage);

        let (alpha, bias) = L1[10].clone();
        let (x, y, z) = line_evaluation(alpha, bias, P2_sage);
        assert_eq!(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime(
                    "19988749844081165334827762793814350502116891030114757846662472163944853248859"
                )
                    .unwrap(),
                bn256::fq::Fq::from_str_vartime(
                    "16234116220677494193936075252196245954744370864079930711074281260312731718107"
                )
                    .unwrap()
            ),
            x
        );
        assert_eq!(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime(
                    "5849570163562045648115193489845684927548545595151466539881169599010675770645"
                )
                    .unwrap(),
                bn256::fq::Fq::from_str_vartime(
                    "3130285240608014807709889950683859105844651234758755309521365436451889463056"
                )
                    .unwrap()
            ),
            y
        );
        assert_eq!(
            Fq2::new(
                bn256::fq::Fq::from_str_vartime("0").unwrap(),
                bn256::fq::Fq::from_str_vartime(
                    "3257012016312567276892866616568971437183359845703546812812284257267928868651"
                )
                    .unwrap()
            ),
            z
        );

        let actual = multi_miller_loop_sage(&[(&P1_sage, L1), (&P2_sage, L2)]);
        assert_eq!(Fq12::new(
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("2579125451910549598126799194972337195136172508947177052057405794715481654051").unwrap(),
                    bn256::fq::Fq::from_str_vartime("8799627207500576818198646407402500173938948809107278073546141842223550815443").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("21254711635776122729543082146190615551870716159044063308362385251619490821191").unwrap(),
                    bn256::fq::Fq::from_str_vartime("15172894816579948629874409203718892278079979405950693229989141494477041252567").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("20119155179103325631029860919524403363393719082193349772072505547385450716052").unwrap(),
                    bn256::fq::Fq::from_str_vartime("16463560053605328069125708406506825705777480445393953022661996923865546566275").unwrap()
                )
            ),
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("11737433663001962884210355938390957295663409438189206545904631195229647825404").unwrap(),
                    bn256::fq::Fq::from_str_vartime("2352394207287859308016976440085523686637738296499997552170469196818269460903").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("11758085804603788208432086676910821478101610524034131161364090902621687492913").unwrap(),
                    bn256::fq::Fq::from_str_vartime("9946296189927120863548066442000400550399885763102241167443723181990702343901").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("20016452407421048126171105605380529407885824070164997819740419654108854409334").unwrap(),
                    bn256::fq::Fq::from_str_vartime("18967240658201571614605289533074013128288058944608046715336000780077678673701").unwrap()
                )
            )
        ), actual);
    }

    #[test]
    fn test_multi_miller_loop_sage() {
        let P1 = g1_scalar_mul(g1_sage(), BigInt::from(87i32));
        let P2 = g1_scalar_mul(g1_sage(), BigInt::from(134i32));

        let Q1 = g2_scalar_mul(g2_sage(), BigInt::from(134i32));
        let Q2 = g2_negate(g2_scalar_mul(g2_sage(), BigInt::from(87i32)));

        let L1 = line_function(Q1);
        let L2 = line_function(Q2);

        let actual = multi_miller_loop_sage(&[(&P1, L1), (&P2, L2)]);
        assert_eq!(Fq12::new(
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("14979159024391170005639198452819736466654857283214903755083678356758517136592").unwrap(),
                    bn256::fq::Fq::from_str_vartime("19179255212017640626421950110077442332223637319173446355334494029062723511588").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("5153843813692294978752855649018278829839715141911891306173539636896001583467").unwrap(),
                    bn256::fq::Fq::from_str_vartime("882594320831470586023090324556522680111944454995811988621287974399994365236").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("6375494376791364655597588792964347300774889947372301757612904362154190445252").unwrap(),
                    bn256::fq::Fq::from_str_vartime("7411177602546289276998323564732516132202756502639327046081197595425717358315").unwrap()
                )
            ),
            Fq6::new(
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("18537294362773696139452051147056860344601322877763095254369218889185173619948").unwrap(),
                    bn256::fq::Fq::from_str_vartime("6501301362025914733678307551845740401477962143257795738951827205743710824435").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("2638916684912325257814418923235928592796227957144644709099582618603445099971").unwrap(),
                    bn256::fq::Fq::from_str_vartime("8497882020477963042360273164135335484186094806658504072297719648400981694791").unwrap()
                ),
                Fq2::new(
                    bn256::fq::Fq::from_str_vartime("20330300560423193451197480003744561686727292889019608644722757092740994444742").unwrap(),
                    bn256::fq::Fq::from_str_vartime("16359697216046849917011637578628027695662599966346456216280707466374786375606").unwrap()
                )
            )
        ), actual);
    }

    fn fq6_negative_of(a: Fq6) -> Fq6 {
        let mut out = a.clone();
        out.c0 = fq2_negative_of(a.c0);
        out.c1 = fq2_negative_of(a.c1);
        out.c2 = fq2_negative_of(a.c2);
        out
    }

    fn fq6_inverse(a: Fq6) -> Fq6 {
        // Algorithm 17

        let XX = fq2_mul(a.c0, a.c0);
        let YY = fq2_mul(a.c1, a.c1);
        let ZZ = fq2_mul(a.c2, a.c2);

        let XY = fq2_mul(a.c0, a.c1);
        let XZ = fq2_mul(a.c0, a.c2);
        let YZ = fq2_mul(a.c1, a.c2);

        let A = ZZ - fq2_mul_beta(XY);
        let B = fq2_mul_beta(XX) - YZ;
        let C = YY - XZ;

        let mut F = fq2_mul_beta(fq2_mul(C, a.c1));
        F += fq2_mul(A, a.c2);
        F += fq2_mul_beta(fq2_mul(B, a.c0));
        F = fq2_inverse(F);

        let c_x = fq2_mul(C, F);
        let c_y = fq2_mul(B, F);
        let c_z = fq2_mul(A, F);

        Fq6::new(c_x, c_y, c_z)
    }

    fn fq12_inverse(a: Fq12) -> Fq12 {
        let mut out = a.clone();
        out.c0 = fq6_negative_of(a.c0);

        let mut t1 = fq6_mul(a.c0, a.c0);
        let mut t2 = fq6_mul(a.c1, a.c1);
        t1 = fq6_mul_tau(t1);
        t1 = t2 - t1;
        t2 = fq6_inverse(t1);

        out.c0 = fq6_mul(out.c0, t2);
        out.c1 = fq6_mul(out.c1, t2);

        out
    }

    // TODO: check if multiplications can be replaced by more efficient square, double, etc.
    fn fq12_exp(a: Fq12, exp: BigInt) -> Fq12 {
        let mut e = to_naf(exp);
        e.reverse();
        e = e[1..].to_vec();

        let mut R = a.clone();
        for (i, kb) in e.into_iter().enumerate() {
            R = fq12_mul(R, R); // square

            if kb == 1 {
                R = fq12_mul(R, a);
            } else if kb == -1 {
                R = fq12_mul(R, fq12_inverse(a));
            }
        }
        R
    }

    fn rth_root(r: BigInt, r_co: BigInt, f: Fq12) -> Fq12 {
        assert_eq!(fq12_exp(f, r_co.clone()), fq12_swap(Fq12::one()));
        let g = r.extended_gcd(&r_co);

        assert_eq!(g.gcd, BigInt::one());
        let r_inv: BigInt = if g.x.is_negative() {
            r_co.clone() + g.x
        } else {
            g.x
        };
        assert_eq!(r.clone().mul(r_inv.clone()).mod_floor(&r_co), BigInt::one());

        let root = fq12_exp(f, r_inv);
        assert_eq!(fq12_exp(root, r), f);

        root
    }

    // https://eprint.iacr.org/2009/457.pdf
    fn tonelli_shanks3(c: Fq12, s: u32, t: BigInt, a: Fq12, k: BigInt) -> Fq12 {
        let mut r = fq12_exp(a, t.clone());

        // compute cubic root of (a^t)^-1, h
        let mut h = fq12_swap(Fq12::one());
        let cc = fq12_exp(c, BigInt::from(3).pow(s - 1));
        let mut c = fq12_inverse(c);

        let mut d = Fq12::zero();
        for i in 1..s {
            let delta = s - i - 1;
            let d = fq12_exp(r, BigInt::from(3).pow(delta));
            if d == cc {
                h = fq12_mul(h, c);
                r = fq12_mul(r, fq12_mul(fq12_mul(c, c), c))
            } else if d == fq12_mul(cc, cc) {
                h = fq12_mul(h, fq12_mul(c, c));
                let c3 = fq12_mul(fq12_mul(c, c), c);
                r = fq12_mul(r, fq12_mul(c3, c3));
            }

            c = fq12_mul(fq12_mul(c, c), c);
        }

        r = fq12_exp(a, k.clone());
        r = fq12_mul(r, h);

        if t == BigInt::from(3).mul(k) + BigInt::one() {
            r = fq12_inverse(r);
        }

        assert_eq!(fq12_exp(r, BigInt::from(3)), a);
        r
    }

    fn compute_final_exp_witness(f: Fq12) -> (Fq12, Fq12){
        let m = mx_x();
        let h = hx_x();
        let p = px_x();
        let r = rx_x();

        let d = BigInt::gcd(&m, &h);

        let mm = m.div_mod_floor(&d).0;

        assert_eq!(d.clone().mul(mm.clone()).mul(r.clone()), lambdax_x());

        let s = 3; // constant

        let p_pow_12_minus_1 = p.pow(12).sub(&BigInt::one());
        let t = p_pow_12_minus_1.div_mod_floor(&(BigInt::from(3).pow(s))).0;
        let var_3s_minus_1_mul_t = BigInt::from(3).pow(s - 1).mul(t.clone());

        let k = (t.clone().add(BigInt::one())).div_mod_floor(&BigInt::from(3)).0;

        assert_eq!(r.clone().mul(h.clone()), p_pow_12_minus_1);
        assert_eq!(m.clone().mul(r.clone()), lambdax_x());
        assert_eq!(d.clone().mul(mm.clone()), m);
        assert_eq!(d, BigInt::from(3));
        assert_ne!(fq12_exp(f, var_3s_minus_1_mul_t.clone()), fq12_swap(Fq12::one()));

        let mut rng = OsRng;
        let mut z = Fq12::zero();
        let mut legendre = fq12_swap(Fq12::one());

        // looking for a 'legendre' which is not Fq12::one()
        loop {
            z = Fq12::random(&mut rng);
            legendre = fq12_exp(z, var_3s_minus_1_mul_t.clone());
            if fq12_swap(Fq12::one()) != legendre {
                break
            }
        }

        let w = fq12_exp(z, t.clone());
        assert_ne!(w, fq12_swap(Fq12::one()));
        assert_ne!(fq12_exp(w, var_3s_minus_1_mul_t.clone()), fq12_swap(Fq12::one()));
        assert_eq!(fq12_exp(w, h.clone()), fq12_swap(Fq12::one()));

        let mut wi = w;
        if fq12_exp(fq12_mul(f, w), var_3s_minus_1_mul_t.clone()) != fq12_swap(Fq12::one()) {
            assert_eq!(fq12_exp(fq12_mul(f, fq12_mul(w, w)), var_3s_minus_1_mul_t.clone()), fq12_swap(Fq12::one()));
            wi = fq12_mul(w, w)
        }

        let f1 = fq12_mul(f, wi);
        assert_eq!(fq12_exp(f1, h.clone()), fq12_swap(Fq12::one()));

        let f2 = rth_root(r.clone(), h.clone(), f1);
        assert_eq!(fq12_exp(f2, r.clone()), f1);

        let f3 = rth_root(mm.clone(), r.clone().mul( h.clone()), f2);
        assert_eq!(fq12_exp(f3, mm.clone().mul(r.clone())), f1);

        let c= tonelli_shanks3(w, s, t, f3, k);

        (c, wi)
    }

    #[test]
    fn test_pairing_verification_on_sage_input() {
        // sage input
        let P1 = g1_scalar_mul(g1_sage(), BigInt::from(176i32));
        let P2 = g1_scalar_mul(g1_sage(), BigInt::from(19292i32));

        let Q1 = g2_scalar_mul(g2_sage(), BigInt::from(19292i32));
        let Q2 = g2_negate(g2_scalar_mul(g2_sage(), BigInt::from(176i32)));

        // canonical pairing verification
        let mut mml_res = multi_miller_loop(&[(&P1, &g2_swap(Q1)), (&P2, &g2_swap(Q2))]);
        assert_eq!(Gt::identity(), mml_res.final_exponentiation());


        let L1 = line_function(Q1);
        let L2 = line_function(Q2);

        // optimized pairing verification
        let f = multi_miller_loop_sage(&[(&P1, L1), (&P2, L2)]);
        let (c, wi) = compute_final_exp_witness(f);
        assert_eq!(fq12_exp(c, lambdax_x()), fq12_mul(f, wi));
    }

    #[test]
    fn test_pairing_verification_on_canonical_input() {
        // canonical input
        let P1 =  G1Affine::generator();
        let P2 = -G1Affine::generator();

        let Q1 = G2Affine::generator();
        let Q2 = G2Affine::generator();

        // canonical pairing verification
        let mut mml_res = multi_miller_loop(&[(&P1, &Q1), (&P2, &Q2)]);
        assert_eq!(Gt::identity(), mml_res.final_exponentiation());


        let L1 = line_function(g2_swap(Q1));
        let L2 = line_function(g2_swap(Q2));

        // optimized pairing verification
        let f = multi_miller_loop_sage(&[(&P1, L1), (&P2, L2)]);
        let (c, wi) = compute_final_exp_witness(f);
        assert_eq!(fq12_exp(c, lambdax_x()), fq12_mul(f, wi));
    }
}
