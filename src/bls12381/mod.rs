mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;
mod g1;
mod g2;

pub use engine::*;
pub use fq::*;
pub use fq12::*;
pub use fq2::*;
pub use fq6::*;
pub use fr::*;
pub use g1::*;
pub use g2::*;

const BLS_X: [u8; 64] = [
    1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

const ENDO_PARAMS: EndoParameters = EndoParameters {
    // round(b2/n)
    gamma2: [0x63f6e522f6cfee30u64, 0x7c6becf1e01faadd, 0x01, 0x0],
    // round(-b1/n)
    gamma1: [0x02u64, 0x0, 0x0, 0x0],
    b1: [0x01u64, 0x0, 0x0, 0x0],
    b2: [0x0000000100000000, 0xac45a4010001a402, 0x0, 0x0],
};

use crate::arithmetic::{mul_512, sbb, CurveEndo, EndoParameters};
use ff::{PrimeField, WithSmallOrderMulGroup};
crate::endo!(G1, Fr, ENDO_PARAMS);
crate::endo!(G2, Fr, ENDO_PARAMS);

#[cfg(test)]
mod tests {
    use std::ops::{Add, Div, Mul, Neg, Sub, SubAssign};
    use ff::{Field, PrimeField};
    use group::{Curve, Group};
    use group::prime::PrimeCurveAffine;
    use pairing::MillerLoopResult;
    use crate::bls12381::{BLS_X, Fq12, Fq2, Fq6, FROBENIUS_COEFF_FQ12_C1, FROBENIUS_COEFF_FQ6_C2, G1Affine, G2, G2Affine, Gt, multi_miller_loop};
    use num_bigint::{BigInt, BigUint, ToBigInt};
    use num_integer::{div_floor, Integer};
    use num_traits::real::Real;
    use num_traits::{FromPrimitive, Num, One, Signed, ToPrimitive, Zero};
    use pasta_curves::arithmetic::CurveAffine;
    use pasta_curves::Fq;
    use rand_core::OsRng;
    use crate::bls12381;
    use crate::ff_ext::ExtField;

    #[test]
    fn test_canonical_pairing_example() {
        let a_g1 = G1Affine::generator().mul(bls12381::fr::Fr::from(10)).to_affine();
        let a_g2 = G2Affine::generator().mul(bls12381::fr::Fr::from(20)).to_affine();

        let b_g1 = -(G1Affine::generator().mul(bls12381::fr::Fr::from(20)).to_affine());
        let b_g2 = G2Affine::generator().mul(bls12381::fr::Fr::from(10)).to_affine();

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



    // curve / optimisation parameters
    fn x() -> BigInt {
        BigInt::from_i128(-0xd201000000010000).unwrap()
    }

    fn modulus() -> BigInt {
        BigInt::from_str_radix(
            "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
            16,
        )
            .unwrap()
    }

    fn modulus_minus_one() -> BigInt {
        BigInt::from_str_radix(
            "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa",
            16,
        )
            .unwrap()
    }

    fn qx_x() -> BigInt {
        BigInt::from_str_radix(
            "793479390729215512464072076108042875612148084278281638540369054903177558538675143759988661420032",
            10,
        )
            .unwrap()
    }

    fn rx_x() -> BigInt {
        BigInt::from_str_radix(
            "52435875175126190479447740508185965837690552500527637822603658699938581184513",
            10,
        )
            .unwrap()
    }

    fn hx_x() -> BigInt {
        BigInt::from_str_radix(
            "1187953094646949460964978789086633769667870913708383532467360559333596611577789496090610146346167525077031903026387460297602145765838579871154117300455841852190629055205186425600590251150517515367911443174768705320573932810930198120802599044705717213480718499961914049370800645289059459984836450604518002989477047170323071258342970036221874890129240095699552337979714210019009615509730638623000191466067385755034978291908382510289449381875044778699119385930236387772922692338048126021229443474326936036965826396020945253945759415758850617310018220324301384943871562440814222536907575606067687182111423556285515824827520364666967830224662633247787077074782720728131656772547794195854262675414567816583080110720835578359169720710074884593917028861514419257932964751212724467326773344978244858101173043370146651784103401701305804523361399056259184066219923820862438873422105479651196757069698007559301663313367301644956596508010228054120661886621909379319031507296096523254415667360258340428512229113253229789201650038514272520862731499154749826274492214922606589068396856344575",
            10,
        )
            .unwrap()
    }

    fn lambdax_x() -> BigInt {
        BigInt::from_str_radix(
            "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129030796414117214202539",
            10,
        )
            .unwrap()
    }

    fn mx_x() -> BigInt {
        BigInt::from_str_radix(
            "76329603384216526031706109802092473003",
            10,
        )
            .unwrap()
    }

    #[test]
    fn test_miller_loop() {
        // canonical input
        let P1 =  G1Affine::generator();
        let P2 = -G1Affine::generator();

        let Q1 = G2Affine::generator();
        let Q2 = G2Affine::generator();

        // canonical pairing verification
        let mut mml_res = multi_miller_loop(&[(&P1, &Q1), (&P2, &Q2)]);
        assert_eq!(Gt::identity(), mml_res.final_exponentiation());

        let f1 = miller_loop(P1, g2_swap(Q1));
        let f2 = miller_loop(P2, g2_swap(Q2));
        let f = fq12_mul(f1, f2);

        let (c, wi) = compute_final_exp_witness(f);
    }


    fn compute_final_exp_witness(f: Fq12) -> (Fq12, Fq12){
        let m = mx_x();
        let h = hx_x();

        // TODO: Why not working? Most likely, either fq12_exp or input f is not correct
        //assert_eq!(fq12_exp(f, h.clone()), fq12_swap(Fq12::one()));

        let q = qx_x();
        let r = rx_x();
        let d = BigInt::gcd(&m, &h);

        let mm = m.div_mod_floor(&d).0;

        assert_eq!(d.clone() * (mm.clone() * r.clone()), lambdax_x());

        // check: r(x) | q(x)^4 − q(x)^2 + 1
        assert!((modulus().pow(4) - modulus().pow(2) + BigInt::one()).is_multiple_of(&r));
        assert!((qx_x().pow(4) - qx_x().pow(2) + BigInt::one()).is_multiple_of(&r));

        assert_eq!(m.clone() * r.clone(), lambdax_x());

        let q_pow_12_minus_1 = q.pow(12).sub(&BigInt::one());
        assert_eq!(r.clone() * h.clone(), q_pow_12_minus_1);

        // check: gcd(λ, q^12 − 1) = d * r = gcd(m, h) * r
        assert_eq!(BigInt::gcd(&lambdax_x(), &q_pow_12_minus_1), BigInt::from(d.clone()) * r.clone());

        assert_eq!(lambdax_x(), r.clone() * mm.clone() * d.clone());

        // ...

        (Fq12::one(), Fq12::one())
    }

    // TODO: check if multiplications can be replaced by more efficient square, double, etc.
    // TODO: find a way to test the correctness
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
    fn fq6_mul_tau(a: Fq6) -> Fq6 {
        let tx = a.c1;
        let ty = a.c2;
        let tz = fq2_mul_beta(a.c0);
        Fq6::new(tx, ty, tz)
    }

    fn fq2_mul_beta(a: Fq2) -> Fq2 {
        // FIXME: probably this function is problematic

        //fq2_swap(Fq2::mul_by_nonresidue(&fq2_swap(a)))
        let tx = a.c0 + a.c1;
        let ty = a.c1 - a.c0;
        Fq2::new(tx, ty)

        /*let out = fq2_swap(a.clone());
        out.mul_by_nonresidue();
        fq2_swap(out)*/
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

    fn fq6_mul(a: Fq6, b: Fq6) -> Fq6 {
        let out = fq6_swap(a) * fq6_swap(b);
        fq6_swap(out)
    }

    fn fq6_negative_of(a: Fq6) -> Fq6 {
        let mut out = a.clone();
        out.c0 = fq2_negative_of(a.c0);
        out.c1 = fq2_negative_of(a.c1);
        out.c2 = fq2_negative_of(a.c2);
        out
    }

    fn miller_loop(p: G1Affine, q: G2Affine) -> Fq12 {
        let mut f = Fq12::new(
            Fq6::new(
                Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::zero()),
                Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::zero()),
                Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::zero()),
            ),
            Fq6::new(
                Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::zero()),
                Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::zero()),
                Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::one()),
            ),
        );

        let mQ = g2_negate(q);
        let mut T = q.clone();
        let mut T_projective = T.to_curve();
        T_projective.z = fq2_swap(T_projective.z);

        let Qp = fq2_mul(q.y, q.y);

        let mut a: Fq2 = Fq2::zero();
        let mut b: Fq2 = Fq2::zero();
        let mut c: Fq2 = Fq2::zero();


        // use negated x as 'to_naf' doesn't support negative values
        let mut naf = to_naf(x().neg());
        naf.reverse();
        naf = naf[1..].to_vec();

        // TODO: Probably we need BLS_X using instead of naf
        for (i, naf_i) in naf.into_iter().enumerate() {
            f = fq12_mul(f, f);

            (a, b, c, T, T_projective) = line_func_double(T_projective, p);

            f = fq12_mul_line_base(f, (a, b, c));

            if naf_i == 1 {
                (a, b, c, T, T_projective) = line_func_add(T_projective, q, p, Qp);
                f = fq12_mul_line_base(f, (a, b, c));
            } else if naf_i == -1 {
                (a, b, c, T, T_projective) = line_func_add(T_projective, mQ, p, Qp);
                f = fq12_mul_line_base(f, (a, b, c));
            }
        }

        assert_eq!(T, g2_scalar_mul(q, BigInt::from(x().neg())));

        // since x is negative, we need to conjugate the output
        fq12_conjugate(f)
    }

    // functionality
    fn fq2_conjugate(a: Fq2) -> Fq2 {
        let mut out = fq2_swap(a.clone());
        out.conjugate();
        fq2_swap(out)
    }

    fn to_naf(x: BigInt) -> Vec<i32> {
        assert!(x.is_positive());

        let mut x = x.clone();

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

        z.iter().map(|item| item.to_i32().unwrap()).collect::<Vec<i32>>()
    }

    fn g2_scalar_mul(g2: G2Affine, value: BigInt) -> G2Affine {
        let mut e = to_naf(value);
        e.reverse();
        e = e[1..].to_vec();

        let g2_initial = g2.clone();
        let mut result = g2_initial.clone();

        for (i, kb) in e.into_iter().enumerate() {
            result = g2_affine_double(result);
            if kb == 1 {
                result = g2_affine_add(result, g2_initial);
            } else if kb == -1 {
                result = g2_affine_add(result, g2_negate(g2_initial));
            }
        }

        result
    }

    fn g2_affine_double(a: G2Affine) -> G2Affine {
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
        let c_z = fq2_mul(a.y, Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::one()));
        let c_z = c_z.double();

        g2_force_affine(c_x, c_y, c_z)
    }
    fn g2_affine_add(a: G2Affine, b: G2Affine) -> G2Affine {
        // http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl

        // TODO: handle infinite point cases

        let a_z = Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::one());
        let b_z = Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::one());

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
            return g2_affine_double(a);
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


    fn fq12_mul_line_base(r: Fq12, abc: (Fq2, Fq2, Fq2)) -> Fq12 {
        let fl = Fq12::new(
            Fq6::new(Fq2::zero(), abc.0, abc.1),
            Fq6::new(Fq2::zero(), Fq2::zero(), abc.2),
        );

        fq12_mul(r, fl)
    }

    fn fq12_mul(a: Fq12, b: Fq12) -> Fq12 {
        let out = fq12_swap(a) * fq12_swap(b);
        fq12_swap(out)
    }

    fn fq2_mul(a: Fq2, b: Fq2) -> Fq2 {
        let out = fq2_swap(a) * fq2_swap(b);
        fq2_swap(out)
    }

    fn g2_negate(a: G2Affine) -> G2Affine {
        let mut out = a.clone();
        out.y = fq2_negative_of(out.y);
        out
    }

    fn fq2_negative_of(a: Fq2) -> Fq2 {
        let mut a_copy = a.clone();
        a_copy.c0 = fq_additive_inverse(a.c0);
        a_copy.c1 = fq_additive_inverse(a.c1);
        a_copy
    }

    fn fq_additive_inverse(a: bls12381::fq::Fq) -> bls12381::fq::Fq {
        let mod_1 =
            bls12381::fq::Fq::from_str_vartime(modulus_minus_one().to_str_radix(10).as_str()).unwrap();
        let mut output = mod_1 - a;
        output += bls12381::fq::Fq::one();
        output
    }

    fn line_func_double(r: G2, q: G1Affine) -> (Fq2, Fq2, Fq2, G2Affine, G2){
        let r_t = fq2_mul(r.z, r.z);
        let A = fq2_mul(r.x, r.x);
        let B = fq2_mul(r.y, r.y);
        let C = fq2_mul(B, B);
        let D = r.x + B;
        let D = fq2_mul(D, D);
        let D = D - A;
        let D = D - C;
        let D = D + D;
        let E = A + A + A;
        let F = fq2_mul(E, E);
        let C8 = C + C + C + C + C + C + C + C;
        let r_x = F - D - D;
        let r_y = fq2_mul(E, D - r_x) - C8;
        let r_z = fq2_mul(r.y + r.z, r.y + r.z) - B - r_t;

        assert_eq!(r_z, fq2_mul(r.y, r.z + r.z));

        let mut r_out = g2_force_affine(r_x, r_y, r_z);
        // we use is_on_curve from halo2curves, so swapping r_out
        assert_eq!(g2_swap(r_out).is_on_curve().unwrap_u8(), 0x01);

        let a = r.x + E;
        let a = fq2_mul(a, a);
        let tmp = A + F + B + B + B + B;
        let a = a - tmp;
        let t = fq2_mul(E, r_t);
        let t = t + t;
        let b = fq2_negative_of(t);
        let b = fq2_mul_scalar(b, q.x);

        let c = fq2_mul(r_z, r_t);
        let c = c + c;
        let c = fq2_mul_scalar(c, q.y);

        let aux_inv = fq2_mul(r_t, r_z);
        let aux_inv = aux_inv + aux_inv;
        let aux_inv = fq2_inverse(aux_inv);

        let a = fq2_mul(a, aux_inv);
        let b = fq2_mul(b, aux_inv);
        let c = fq2_mul(c, aux_inv);

        let mut projective = G2::identity();
        projective.x = r_x;
        projective.y = r_y;
        projective.z = r_z;

        (a, b, c, r_out, projective)

    }

    fn fq2_inverse(a: Fq2) -> Fq2 {
        let t = a.c0.square() + a.c1.square();
        let inv = t.invert().unwrap();
        let c_x = fq_additive_inverse(a.c0) * inv;
        let c_y = a.c1 * inv;
        Fq2::new(c_x, c_y)
    }

    fn fq2_mul_scalar(a: Fq2, scalar: bls12381::fq::Fq) -> Fq2 {
        let mut a_copy = a.clone();
        a_copy.c0 = a.c0 * scalar;
        a_copy.c1 = a.c1 * scalar;
        a_copy
    }

    fn g2_force_affine(x: Fq2, y: Fq2, z: Fq2) -> G2Affine {
        if z == Fq2::zero()
            || z == Fq2::new(bls12381::fq::Fq::zero(), bls12381::fq::Fq::one())
            || z == Fq2::new(bls12381::fq::Fq::one(), bls12381::fq::Fq::zero())
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
    fn line_func_add(r: G2, p: G2Affine, q: G1Affine, r2: Fq2) -> (Fq2, Fq2, Fq2, G2Affine, G2) {
        let r_t = fq2_mul(r.z, r.z);
        let B = fq2_mul(p.x, r_t);
        let D = p.y + r.z;
        let D = fq2_mul(D, D);
        let D = D - r2;
        let D = D - r_t;
        let D = fq2_mul(D, r_t);

        let H = B - r.x;
        let I = fq2_mul(H, H);

        let E = I + I + I + I;
        let J = fq2_mul(H, E);
        let L1 = D - r.y;
        let L1 = L1 - r.y;

        let V = fq2_mul(r.x, E);

        let r_x = fq2_mul(L1, L1);
        let r_x = r_x - J;
        let r_x = r_x - V - V;

        let r_z = r.z + H;
        let r_z = fq2_mul(r_z, r_z);
        let r_z = r_z - r_t;
        let r_z = r_z - I;

        let t = V - r_x;
        let t = fq2_mul(t, L1);
        let t2 = fq2_mul(r.y, J);
        let t2 = t2 + t2;
        let r_y = t - t2;

        //print_fq2(r_x);
        //print_fq2(r_y);
        //print_fq2(r_z);
        //println!();

        let r_out = g2_force_affine(r_x, r_y, r_z);

        let t = p.y + r_z;
        let t = fq2_mul(t, t);
        let t = t - r2;
        let t = t - fq2_mul(r_z, r_z);

        let t2 = fq2_mul(L1, p.x);
        let t2 = t2 + t2;
        let a = t2 - t;

        let c = fq2_mul_scalar(r_z, q.y);
        let c = c + c;

        let b = fq2_negative_of(L1);
        let b = fq2_mul_scalar(b, q.x);
        let b = b + b;

        let aux_inv = fq2_inverse(r_z + r_z);
        let a = fq2_mul(a, aux_inv);
        let b = fq2_mul(b, aux_inv);
        let c = fq2_mul(c, aux_inv);

        let mut projective = G2::identity();
        projective.x = r_x;
        projective.y = r_y;
        projective.z = r_z;

        (a, b, c, r_out, projective)
    }
    fn fq2_exp(a: Fq2, exp: BigInt) -> Fq2 {
        let mut e = to_naf(exp);
        e.reverse();
        e = e[1..].to_vec();

        let mut R = a.clone();
        for (i, kb) in e.into_iter().enumerate() {
            R = fq2_mul(R, R); // square

            if kb == 1 {
                R = fq2_mul(R, a);
            } else if kb == -1 {
                R = fq2_mul(R, fq2_inverse(a));
            }
        }
        R
    }
    fn fq12_conjugate(a: Fq12) -> Fq12 {
        let mut out = fq12_swap(a.clone());
        out.conjugate();
        fq12_swap(out)
    }
}
