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
    use crate::bls12381;
    use crate::bls12381::{multi_miller_loop, Fq12, G1Affine, G2Affine, Gt};
    use ff::Field;
    use group::Curve;
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{Num, One, Signed};
    use pairing::MillerLoopResult;
    use rand::Rng;
    use rand_core::OsRng;
    use std::ops::{Div, Mul, Neg, Sub};

    fn qx_x() -> BigInt {
        BigInt::from_str_radix(
            "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787",
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
            "322277361516934140462891564586510139908379969514828494218366688025288661041104682794998680497580008899973249814104447692778988208376779573819485263026159588510513834876303014016798809919343532899164848730280942609956670917565618115867287399623286813270357901731510188149934363360381614501334086825442271920079363289954510565375378443704372994881406797882676971082200626541916413184642520269678897559532260949334760604962086348898118982248842634379637598665468817769075878555493752214492790122785850202957575200176084204422751485957336465472324810982833638490904279282696134323072515220044451592646885410572234451732790590013479358343841220074174848221722017083597872017638514103174122784843925578370430843522959600095676285723737049438346544753168912974976791528535276317256904336520179281145394686565050419250614107803233314658825463117900250701199181529205942363159325765991819433914303908860460720581408201373164047773794825411011922305820065611121544561808414055302212057471395719432072209245600258134364584636810093520285711072578721435517884103526483832733289802426157301542744476740008494780363354305116978805620671467071400711358839553375340724899735460480144599782014906586543813292157922220645089192130209334926661588737007768565838519456601560804957985667880395221049249803753582637708560",
            10,
        )
            .unwrap()
    }

    fn lambdax_x() -> BigInt {
        BigInt::from_str_radix(
            "-970223322989879234235347398217909883513800123614040636629605779369206487123988021762754806436420539225018110660430574629597272117950833190645132350905828054049340852793006395264263886264819612507783037924218149124701826056046922196945843543463363290689079493037063497864477125563880259435400602850050818969493853842620644220539444087494836928552324610182844026424",
            10,
        )
            .unwrap()
    }

    fn mx_x() -> BigInt {
        BigInt::from_str_radix(
            "-18503044332711365414904791911389002046454456235902880068437438477414596808076536185222086137936701081884770111977142993511392347341594134320909625744350654635267487273453992200083136256093079477666275889079364852447945068536815482946199946320650919671611689786576455244427459274102574648",
            10,
        )
            .unwrap()
    }

    const BLS12_381_EMBEDDING_DEGREE: u32 = 12u32;

    fn sample_f() -> Fq12 {
        let mut scalar_1 = bls12381::fr::Fr::from(10);
        let mut scalar_2 = bls12381::fr::Fr::from(20);

        let mut p1 = G1Affine::generator().mul(scalar_1).to_affine();
        let mut q1 = G2Affine::generator().mul(scalar_2).to_affine();
        let mut p2 = -(G1Affine::generator().mul(scalar_2).to_affine());
        let mut q2 = G2Affine::generator().mul(scalar_1).to_affine();

        // canonical pairing verification test
        let mut f = multi_miller_loop(&[(&p1, &q1), (&p2, &q2)]);
        assert_eq!(Gt::identity(), f.final_exponentiation());

        // BLS12-381 curve parameters according to https://eprint.iacr.org/2024/640.pdf (Section 4.3.1)
        let mut rng = OsRng;
        let m = mx_x();
        let h = hx_x();
        let d = BigInt::gcd(&m, &h);
        let k = BLS12_381_EMBEDDING_DEGREE;
        let q_pow_12_minus_1 = qx_x().pow(k).sub(&BigInt::one());

        // According to https://eprint.iacr.org/2024/640.pdf (Section 4.2.2), the output of Miller loop doesn't
        // require scaling if it is full λ-th residue, which means that it is λ/d-th residue (which is always true)
        // and d-residue (which is not always true). So we re-generate new Miller loop input until we get f
        // which is d-residue. This means that scale factor is 1, so we don't need to implement the scaling algorithm
        // now for simplicity.
        while Fq12::one() != fq12_exp_halo2(f, q_pow_12_minus_1.div_floor(&d)) {
            scalar_1 = bls12381::fr::Fr::from(rng.gen::<u64>());
            scalar_2 = bls12381::fr::Fr::from(rng.gen::<u64>());
            p1 = G1Affine::generator().mul(scalar_1).to_affine();
            q1 = G2Affine::generator().mul(scalar_2).to_affine();
            p2 = -(G1Affine::generator().mul(scalar_2).to_affine());
            q2 = G2Affine::generator().mul(scalar_1).to_affine();
            let mml_res = multi_miller_loop(&[(&p1, &q1), (&p2, &q2)]);
            assert_eq!(Gt::identity(), mml_res.final_exponentiation());
            f = mml_res;
        }

        f
    }

    #[test]
    fn test_miller_loop_using_python_script() {
        let f = sample_f();

        let h3 = BigInt::from_str_radix("2366356426548243601069753987687709088104621721678962410379583120840019275952471579477684846670499039076873213559162845121989217658133790336552276567078487633052653005423051750848782286407340332979263075575489766963251914185767058009683318020965829271737924625612375201545022326908440428522712877494557944965298566001441468676802477524234094954960009227631543471415676620753242466901942121887152806837594306028649150255258504417829961387165043999299071444887652375514277477719817175923289019181393803729926249507024121957184340179467502106891835144220611408665090353102353194448552304429530104218473070114105759487413726485729058069746063140422361472585604626055492939586602274983146215294625774144156395553405525711143696689756441298365274341189385646499074862712688473936093315628166094221735056483459332831845007196600723053356837526749543765815988577005929923802636375670820616189737737304893769679803809426304143627363860243558537831172903494450556755190448279875942974830469855835666815454271389438587399739607656399812689280234103023464545891697941661992848552456326290792224091557256350095392859243101357349751064730561345062266850238821755009430903520645523345000326783803935359711318798844368754833295302563158150573540616830138810935344206231367357992991289265295323280", 10).unwrap();
        let s = BigInt::from_str_radix("13871344871029839191", 10).unwrap();
        let e = BigInt::from_str_radix("1831641335620623066030493719814750505730353469232607408984136636587602798828541096315617787337660580896596326998534644457004429745922908233143891672764988414950129738729858091158205555565902620547217447631024516567208419316911441452506166124582175761417927641295939217754840037829611255709087158279242957940625143296408819749143295856831426628369017685600998542951914342403397019136966636888082401920580506454124412301032882712464981559302655280563177865634616387141297622595244784621594534214052542425250197179149807006335538265574872503469097186776521644002111433399650909641885318866460936246236243529502790661531948114436384119347931262823673992912046778604767497085407555474313407712735307333317882930639512124999103536572343436126871997578016786488693723339409687453114134784227807901874533332546112824417573327256689835936839706774393384112352709458439683933327013619962098977697684860933918339286553214389828566701676847442687594385334748280300467909322727863382808795828794016038705718095323209413819963715797469895587732751744906760424173767179546312457734136109353976884272162867417651164377486191418315758965239169722140695318540646880198727657740603762882647173358897627586506579349681824692808439213054305286657596509971558805903962977863660242510377179301897268259", 10).unwrap();

        let x = fq12_exp_halo2(f, h3);
        let shift = fq12_exp_halo2(x, s);
        let c = fq12_exp_halo2(shift * f, e);


        // FIXME: Why not working (prints FALSE)?
        println!("{}", fq12_exp_halo2(c, lambdax_x()) == f * shift);

        let mut lambda: [u64; 8] = [
            0x0000000000000400,
            0x2409555221667393,
            0x4177898257359041,
            0x5655688281993900,
            0x7885332058136124,
            0x0316504908378644,
            0x4268762912903079,
            0x6414117214202539,
        ];

        println!("{}", c.pow_vartime(lambda.to_vec()) == f * shift);

        lambda.reverse();

        println!("{}", c.pow_vartime(lambda.to_vec()) == f * shift);
    }

    #[test]
    fn test_miller_loop() {
        let f = sample_f();
        compute_final_exp_witness(f);
    }

    fn fq12_exp_halo2(f: Fq12, e: BigInt) -> Fq12 {
        let exp = e.clone().to_u64_digits();
        let mut power = f.pow_vartime(exp.1);
        if e.is_negative() {
            power.conjugate()
        }
        power
    }

    fn compute_final_exp_witness(f: Fq12) -> (Fq12, Fq12) {
        // BLS12-381 curve parameters according to https://eprint.iacr.org/2024/640.pdf (Section 4.3.1)
        let m = mx_x();
        let h = hx_x();
        let q = qx_x();
        let r = rx_x();
        let d = BigInt::gcd(&m, &h);
        assert_eq!(d, BigInt::from(8));
        let m_mark = m.div_floor(&d);

        // Assertions useful for debugging and testing validity of BLS12-381 parameters computed,
        // according to https://eprint.iacr.org/2024/640.pdf (Section 4.3):
        //
        // 1) λ = d * r * m'
        // 2) r(x) | q(x)^4 − q(x)^2 + 1
        // 3) m * r = λ
        // 4) r * h = q^12 - 1
        // 5) gcd(λ, q^12 − 1) = d * r

        assert_eq!(d.clone() * (m_mark.clone() * r.clone()), lambdax_x());
        assert!((qx_x().pow(4) - qx_x().pow(2) + BigInt::one()).is_multiple_of(&r));
        assert_eq!(m.clone() * r.clone(), lambdax_x());
        let q_pow_12_minus_1 = q.pow(12).sub(&BigInt::one());
        assert_eq!(r.clone() * h.clone(), q_pow_12_minus_1);
        assert_eq!(
            BigInt::gcd(&lambdax_x(), &q_pow_12_minus_1),
            BigInt::from(d.clone()) * r.clone()
        );

        // check that no scaling is required: e.g. f is d-th residue
        assert_eq!(
            Fq12::one(),
            fq12_exp_halo2(f, q_pow_12_minus_1.div_floor(&d))
        );

        // Computing r-th root, according to https://eprint.iacr.org/2024/640.pdf (Section 4.3.2)
        let r_mark = r.clone().modinv(&h.clone()).unwrap();
        let c = fq12_exp_halo2(f, r_mark.clone());

        // Assertions for debugging:
        //
        // 1) r' * r = 1 mod h
        // 2) f^h = 1
        // 3) f^h = c^(r * h)
        // 4) c^r = f
        // 5) c^r = f^(r * r')
        //
        // 6) c^r = f^(1 + h * s) !!! (works with s == d)
        // 7) r * r′= 1 + h * s !!! (doesn't work)

        assert_eq!(
            (r_mark.clone() * r.clone()).mod_floor(&h.clone()),
            BigInt::one()
        );
        assert_eq!(fq12_exp_halo2(f, h.clone()), Fq12::one());
        assert_eq!(
            fq12_exp_halo2(f, h.clone()),
            fq12_exp_halo2(c, r.clone() * h.clone())
        );
        assert_eq!(fq12_exp_halo2(c, r.clone()), f);
        assert_eq!(
            fq12_exp_halo2(c, r.clone()),
            fq12_exp_halo2(f, r_mark.clone() * r.clone())
        );

        // FIXME: is it correct that s == d?
        let s = d.clone();
        assert_eq!(
            fq12_exp_halo2(c, r.clone()),
            fq12_exp_halo2(f, BigInt::one() + h.clone() * s.clone())
        );

        // FIXME: Why rr′= 1 + h * s not working?
        //assert_eq!(r.clone() * r_mark.clone(), BigInt::one() + h.clone() * s.clone());

        // Computing m'-th root, according to https://eprint.iacr.org/2024/640.pdf (Section 4.3.2)
        // using r-th root as a base while exponentiation, according to algorithm 5 (Section 4.3.2)
        let f1 = c;
        let m_mark_mark = m_mark.clone().modinv(&q_pow_12_minus_1.clone()).unwrap();
        let c = fq12_exp_halo2(f1, m_mark_mark.clone());

        // Assertions for debugging:
        //
        // 1) m'' * m' = 1 mod (q^12 - 1)
        // 2) gcd(m', q^12 - 1) = 1
        // 3) c^(m' * (q^12 - 1)) = 1
        // 4) f1^(q^12 - 1) = 1
        // 5) f1^(q^12 - 1) = c^(m' * (q^12 - 1))
        // 6) c^m' = f1^(m'' * m')
        //
        // 7) c^m' = f1 !!! (doesn't work)
        // 8) c^m' = f1^(1 + (q^12 - 1) * s) (doesn't work with s = d or any other)

        assert_eq!(
            (m_mark_mark.clone() * m_mark.clone()).mod_floor(&q_pow_12_minus_1.clone()),
            BigInt::one()
        );
        assert_eq!(BigInt::gcd(&m_mark, &q_pow_12_minus_1), BigInt::one());
        assert_eq!(
            fq12_exp_halo2(c, m_mark.clone() * q_pow_12_minus_1.clone()),
            Fq12::one()
        );
        assert_eq!(fq12_exp_halo2(f1, q_pow_12_minus_1.clone()), Fq12::one());
        assert_eq!(
            fq12_exp_halo2(f1, q_pow_12_minus_1.clone()),
            fq12_exp_halo2(c, m_mark.clone() * q_pow_12_minus_1.clone())
        );
        assert_eq!(
            fq12_exp_halo2(c, m_mark.clone()),
            fq12_exp_halo2(f1, m_mark_mark.clone() * m_mark.clone())
        );
        // TODO: Why c^m' = f1 not working?
        //assert_eq!(fq12_exp_halo2(c, m_mark.clone()), f1);
        // TODO: Why c^m' = f1^(1 + s*(q^12 - 1)) not working?
        //assert_eq!(fq12_exp_halo2(c, m_mark.clone()), fq12_exp_halo2(f1, BigInt::one() + q_pow_12_minus_1.clone() * s.clone()));

        let f2 = c;

        // FIXME: Is this thinking correct?:
        // since s = d = 8, we need to find 8-th root of f2 as a final computation, using standard
        // Tonelli-Shanks algorithm (for square root computing) over Fq12 3 times: 2^3 = 8.
        // The valid solution could be any of [c1 .. c8] considering that Tonelli-Shanks may return
        // positive / negative root.

        let (c1, c2) = tonelli_shanks2(f2, q_pow_12_minus_1.clone());

        let (c1_1, c2_1) = tonelli_shanks2(c1, q_pow_12_minus_1.clone());
        let (c1_2, c2_2) = tonelli_shanks2(c2, q_pow_12_minus_1.clone());

        let (c1, c2) = tonelli_shanks2(c1_1, q_pow_12_minus_1.clone());
        let (c3, c4) = tonelli_shanks2(c2_1, q_pow_12_minus_1.clone());
        let (c5, c6) = tonelli_shanks2(c1_2, q_pow_12_minus_1.clone());
        let (c7, c8) = tonelli_shanks2(c2_2, q_pow_12_minus_1.clone());

        // TODO: Why final verification (either of) not working (prints FALSE)?
        println!("{}", fq12_exp_halo2(c1, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c2, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c3, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c4, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c5, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c6, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c7, lambdax_x()) == f);
        println!("{}", fq12_exp_halo2(c8, lambdax_x()) == f);

        (Fq12::one(), Fq12::one())
    }

    fn tonelli_shanks2(a: Fq12, modulus: BigInt) -> (Fq12, Fq12) {
        // https://eprint.iacr.org/2009/457.pdf (Table 1)
        let s = 4u32;
        let t = BigInt::from_str_radix("1056180968766935973298317094727349129720276579917759510275632105953966720143584525467678375534976450538765578465124706401179520414005839512675819430566757069238246610141928465766950637165246500340114927381097983658238299241506380945855858205692180675240935100816412276887355919511779477081170611513708805058441133552560439325619778446871960565396621503861301658091113800948928640284214242475992914085753575857173055662621839773156530378803016600478220868219124280395971765031792093192034346711969404552456862943964618691739592020169614661500743338800537928425866522871209523827993287570376048195133453780027469581769674704736285791356089960432515245354200840691228068050847449004047484685145241930846449004327891865448335365551470768305180478279736432533119250858689402048052010528808419161488243369866004172635940073746421410996722360471595721181110626847792454425776378795618273155108264066912480754174632857155444789697344001430818178601806078975806786901827299497480199226460074625472749789730047847818370228473012329257956625102013653022125849710771535456558600132374386894513930623761887534416490401848696972086625571934897737091488565372964455950445940020894473627263940758515198439855038240783863443471324711765496932853921198772467630672273414262838765009111638435703331753341465537647107023369703426547471614433068386633511218403943777627021345503526960699054970705", 10).unwrap();
        assert!(t.is_odd());
        assert_eq!(modulus, t.clone() * BigInt::from(2u64).pow(s));

        // 15469049 is a quadratic non-residue over bls12-381 Fq12
        let mut b = Fq12::zero();
        for _ in 0..15469049u64 {
            b += Fq12::one();
        }

        let c = fq12_exp_halo2(b, t.clone());
        let mut r = fq12_exp_halo2(a, t.clone());
        let mut h = Fq12::one();

        let mut c = c.invert().unwrap();
        for i in 1..s - 1 {
            let d = fq12_exp_halo2(r, BigInt::from(2u64).pow(s - i - 1));
            if d != Fq12::one() {
                h = h * c;
                r = r * c * c;
                c = c * c;
            }
        }

        let solution_1 = h * fq12_exp_halo2(a, (t.clone() + BigInt::one()).div(2));
        let solution_2 = solution_1.neg();

        (solution_1, solution_2)
    }
}
