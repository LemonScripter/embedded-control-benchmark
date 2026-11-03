//! Fedelzeti (beagyazott-kompatibilis) megvalositas: Kolumban-burok biztonsagi kapu egy dron fuggoleges tengelyere.
//! Tulajdonsagok: fix meretu tombok (nincs heap / Vec az algoritmusban), nincs rekurzio, korlatos ciklusok, f32, csuszo ablak INKREMENTALIS frissitessel
//! (az elojel-tabla oszlopa = az ablak-helye -> minden ablakban minden mintanak fuggetlen, adat-fuggetlen elojele van, a tetel feltetele teljesul).
//! Az algoritmus altalanos T: Num tipusra irt: f32 (idomeres) es Cn (muveletszamlalo, minden + - * / sqrt osszehasonlitas szamolt).
use std::fs;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::time::Instant;

const M1: usize = 99; // m - 1
const NMAX: usize = 40;
const Q: usize = 5;
const GDT: f32 = 0.4905;
const VLIM: f32 = 1.0;
const A_NOM: f32 = 0.97;
const B_NOM: f32 = 0.4905;
const PD_EPS: f32 = 1e-6;

trait Num: Copy + PartialOrd + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self> {
    fn c(x: f32) -> Self;
    fn sqrt(self) -> Self;
    fn f(self) -> f32;
}
impl Num for f32 {
    fn c(x: f32) -> f32 { x }
    fn sqrt(self) -> f32 { f32::sqrt(self) }
    fn f(self) -> f32 { self }
}

static ADD: AtomicU64 = AtomicU64::new(0);
static MUL: AtomicU64 = AtomicU64::new(0);
static DIV: AtomicU64 = AtomicU64::new(0);
static SQR: AtomicU64 = AtomicU64::new(0);
static CMP: AtomicU64 = AtomicU64::new(0);
fn counts() -> [u64; 5] { [ADD.load(Relaxed), MUL.load(Relaxed), DIV.load(Relaxed), SQR.load(Relaxed), CMP.load(Relaxed)] }

#[derive(Copy, Clone)]
struct Cn(f32);
impl Add for Cn { type Output = Cn; fn add(self, o: Cn) -> Cn { ADD.fetch_add(1, Relaxed); Cn(self.0 + o.0) } }
impl Sub for Cn { type Output = Cn; fn sub(self, o: Cn) -> Cn { ADD.fetch_add(1, Relaxed); Cn(self.0 - o.0) } }
impl Mul for Cn { type Output = Cn; fn mul(self, o: Cn) -> Cn { MUL.fetch_add(1, Relaxed); Cn(self.0 * o.0) } }
impl Div for Cn { type Output = Cn; fn div(self, o: Cn) -> Cn { DIV.fetch_add(1, Relaxed); Cn(self.0 / o.0) } }
impl Neg for Cn { type Output = Cn; fn neg(self) -> Cn { ADD.fetch_add(1, Relaxed); Cn(-self.0) } }
impl PartialEq for Cn { fn eq(&self, o: &Cn) -> bool { self.0 == o.0 } }
impl PartialOrd for Cn { fn partial_cmp(&self, o: &Cn) -> Option<std::cmp::Ordering> { CMP.fetch_add(1, Relaxed); self.0.partial_cmp(&o.0) } }
impl Num for Cn {
    fn c(x: f32) -> Cn { Cn(x) }
    fn sqrt(self) -> Cn { SQR.fetch_add(1, Relaxed); Cn(self.0.sqrt()) }
    fn f(self) -> f32 { self.0 }
}

struct Filter<T: Num> {
    n: usize,
    sign: [[i8; NMAX]; M1],
    sv: [T; NMAX], su: [T; NMAX], sy: [T; NMAX],
    slot: usize,
    r: [T; 3], b0: [T; 2], bi: [[T; 2]; M1], ai: [[T; 3]; M1],
    th: [T; 2], w: [T; 3],
    zeta: [[T; 2]; M1], pinv: [[T; 3]; M1], kappa: [T; M1], pd: [bool; M1], n_np: usize,
    ok: bool,
}

impl<T: Num> Filter<T> {
    fn new(n: usize, sign: [[i8; NMAX]; M1]) -> Self {
        let z = T::c(0.0);
        Filter { n, sign, sv: [z; NMAX], su: [z; NMAX], sy: [z; NMAX], slot: 0, r: [z; 3], b0: [z; 2], bi: [[z; 2]; M1], ai: [[z; 3]; M1],
                 th: [z; 2], w: [z; 3], zeta: [[z; 2]; M1], pinv: [[z; 3]; M1], kappa: [z; M1], pd: [false; M1], n_np: M1, ok: false }
    }

    /// uj minta beepitese a csuszo ablakba: O(m) = 99 x 5 szorzas-osszeadas (a regi minta ugyanazt az elojel-oszlopot hasznalta)
    fn push(&mut self, v: T, u: T, y: T) {
        let s = self.slot;
        let (ov, ou, oy) = (self.sv[s], self.su[s], self.sy[s]);
        let (dvv, dvu, duu) = (v * v - ov * ov, v * u - ov * ou, u * u - ou * ou);
        let (dyv, dyu) = (y * v - oy * ov, y * u - oy * ou);
        self.r[0] = self.r[0] + dvv; self.r[1] = self.r[1] + dvu; self.r[2] = self.r[2] + duu;
        self.b0[0] = self.b0[0] + dyv; self.b0[1] = self.b0[1] + dyu;
        for i in 0..M1 {
            let a = if self.sign[i][s] > 0 { T::c(1.0) } else { T::c(-1.0) };
            self.bi[i][0] = self.bi[i][0] + a * dyv; self.bi[i][1] = self.bi[i][1] + a * dyu;
            self.ai[i][0] = self.ai[i][0] + a * dvv; self.ai[i][1] = self.ai[i][1] + a * dvu; self.ai[i][2] = self.ai[i][2] + a * duu;
        }
        self.sv[s] = v; self.su[s] = u; self.sy[s] = y;
        self.slot = (s + 1) % self.n;
    }

    /// theta_LS, W = R^-1/2, es a 99 ellipszoid parameterei (n = 2 zart alak)
    fn derive(&mut self) {
        let invn = T::c(1.0 / self.n as f32);
        let (r11, r12, r22) = (self.r[0] * invn, self.r[1] * invn, self.r[2] * invn);
        let (bv, bu) = (self.b0[0] * invn, self.b0[1] * invn);
        let det = r11 * r22 - r12 * r12;
        if !(det > T::c(1e-9)) { self.ok = false; self.n_np = M1; return; }
        self.th = [(r22 * bv - r12 * bu) / det, (r11 * bu - r12 * bv) / det];
        let s = det.sqrt(); let t = (r11 + r22 + s + s).sqrt();
        let (sq11, sq12, sq22) = ((r11 + s) / t, r12 / t, (r22 + s) / t);
        let dsq = sq11 * sq22 - sq12 * sq12;
        let (w11, w12, w22) = (sq22 / dsq, -sq12 / dsq, sq11 / dsq);
        self.w = [w11, w12, w22];
        let (th0, th1) = (self.th[0], self.th[1]);
        let half = T::c(0.5); let one = T::c(1.0); let eps = T::c(PD_EPS);
        let mut np = 0usize;
        for i in 0..M1 {
            let (a11, a12, a22) = (self.ai[i][0] * invn, self.ai[i][1] * invn, self.ai[i][2] * invn);
            let (e1, e2) = (self.bi[i][0] * invn - (a11 * th0 + a12 * th1), self.bi[i][1] * invn - (a12 * th0 + a22 * th1));
            let (f1, f2) = (w11 * e1 + w12 * e2, w12 * e1 + w22 * e2);
            let (aw11, aw12, aw21, aw22) = (a11 * w11 + a12 * w12, a11 * w12 + a12 * w22, a12 * w11 + a22 * w12, a12 * w12 + a22 * w22);
            let (m11, m12, m22) = (w11 * aw11 + w12 * aw21, w11 * aw12 + w12 * aw22, w12 * aw12 + w22 * aw22);
            let (mm11, mm12, mm22) = (m11 * m11 + m12 * m12, m12 * (m11 + m22), m12 * m12 + m22 * m22);
            let (p11, p12, p22) = (one - mm11, -mm12, one - mm22);
            let hd = (p11 - p22) * half;
            let lam = (p11 + p22) * half - (hd * hd + p12 * p12).sqrt();
            if lam > eps {
                let dp = p11 * p22 - p12 * p12;
                let (pi11, pi12, pi22) = (p22 / dp, -p12 / dp, p11 / dp);
                let (u1, u2) = (m11 * f1 + m12 * f2, m12 * f1 + m22 * f2);
                let (z1, z2) = (-(pi11 * u1 + pi12 * u2), -(pi12 * u1 + pi22 * u2));
                self.zeta[i] = [z1, z2]; self.pinv[i] = [pi11, pi12, pi22];
                self.kappa[i] = f1 * f1 + f2 * f2 - (u1 * z1 + u2 * z2);
                self.pd[i] = true;
            } else { self.pd[i] = false; np += 1; }
        }
        self.n_np = np; self.ok = np < Q;
    }

    /// q-adik legnagyobb tartoszam a z-ter egyseg-iranyaban (cx, cy)
    fn support(&self, cx: T, cy: T) -> T {
        let qe = Q as isize - self.n_np as isize;
        if qe <= 0 { return T::c(f32::INFINITY); }
        let qe = qe as usize;
        let zero = T::c(0.0);
        let mut best = [T::c(f32::NEG_INFINITY); Q];
        for i in 0..M1 {
            if !self.pd[i] { continue; }
            let (z, p) = (self.zeta[i], self.pinv[i]);
            let quad = cx * (p[0] * cx + p[1] * cy) + cy * (p[1] * cx + p[2] * cy);
            let q0 = if quad > zero { quad } else { zero };
            let k0 = if self.kappa[i] > zero { self.kappa[i] } else { zero };
            let h = cx * z[0] + cy * z[1] + (k0 * q0).sqrt();
            if h > best[qe - 1] {
                let mut j = qe - 1;
                while j > 0 && h > best[j - 1] { best[j] = best[j - 1]; j -= 1; }
                best[j] = h;
            }
        }
        best[qe - 1]
    }

    /// a burok tartoszama a theta-ter (pv, pu) iranyaban: max_theta P.(theta - theta_hat)
    fn hx(&self, pv: T, pu: T) -> T {
        let (c1, c2) = (self.w[0] * pv + self.w[1] * pu, self.w[1] * pv + self.w[2] * pu);
        let nrm = (c1 * c1 + c2 * c2).sqrt();
        nrm * self.support(c1 / nrm, c2 / nrm)
    }

    fn f_hi(&self, v: T, u: T) -> T { self.th[0] * v + self.th[1] * u - T::c(GDT) + self.hx(v, u) }
    fn f_lo(&self, v: T, u: T) -> T { -(self.th[0] * v + self.th[1] * u) + T::c(GDT) + self.hx(-v, -u) }

    /// biztonsagi kapu: a legkozelebbi biztonsagos u (|v+| <= V minden hihetonek tartott (a, b)-re). Visszateres: (u, allapot) 0=atenged 1=korlatozott 2=vesz 3=tartalek
    fn gate(&self, v: T, u_nom: T) -> (T, u8) {
        let vl = T::c(VLIM); let zero = T::c(0.0); let three = T::c(3.0);
        let nominal = |v: T, u_nom: T| -> T {
            let (uh, ul) = ((vl + T::c(GDT) - T::c(A_NOM) * v) / T::c(B_NOM), (-vl + T::c(GDT) - T::c(A_NOM) * v) / T::c(B_NOM));
            if u_nom < ul { ul } else if u_nom > uh { uh } else { u_nom }
        };
        if !self.ok { return (nominal(v, u_nom), 3); }
        if self.f_hi(v, u_nom) <= vl && self.f_lo(v, u_nom) <= vl { return (u_nom, 0); }
        let bmin = self.th[1] - self.hx(zero, T::c(-1.0));
        if !(bmin > T::c(0.05)) { return (nominal(v, u_nom), 3); }
        let half = T::c(0.5);
        // u_hi: max u, amelyre f_hi <= V (f_hi u-ban novekvo, ha b > 0 a burokban)
        let u_hi = if self.f_hi(v, zero) > vl { T::c(-1.0) } else if self.f_hi(v, three) <= vl { three } else {
            let (mut lo, mut hi) = (zero, three);
            for _ in 0..8 { let mid = (lo + hi) * half; if self.f_hi(v, mid) <= vl { lo = mid } else { hi = mid } }
            lo };
        let u_lo = if self.f_lo(v, three) > vl { T::c(4.0) } else if self.f_lo(v, zero) <= vl { zero } else {
            let (mut lo, mut hi) = (zero, three);
            for _ in 0..8 { let mid = (lo + hi) * half; if self.f_lo(v, mid) <= vl { hi = mid } else { lo = mid } }
            hi };
        if u_lo <= u_hi { let u = if u_nom < u_lo { u_lo } else if u_nom > u_hi { u_hi } else { u_nom }; return (u, 1); }
        let (mut lo, mut hi) = (zero, three);
        for _ in 0..8 { let mid = (lo + hi) * half; if self.f_hi(v, mid) - self.f_lo(v, mid) <= zero { lo = mid } else { hi = mid } }
        ((lo + hi) * half, 2)
    }
}

struct Tv { n: usize, samples: Vec<[f32; 3]>, sign: [[i8; NMAX]; M1], th: [f32; 2], w: [f32; 3], pd: usize, hx: Vec<[f32; 3]>, gate: Vec<[f32; 3]> }

fn read_tv(path: &str) -> Tv {
    let txt = fs::read_to_string(path).expect("tv fajl");
    let mut it = txt.split_whitespace();
    let mut nx = |it: &mut std::str::SplitWhitespace| it.next().unwrap().to_string();
    assert_eq!(nx(&mut it), "N"); let n: usize = nx(&mut it).parse().unwrap();
    assert_eq!(nx(&mut it), "S");
    let samples: Vec<[f32; 3]> = (0..n).map(|_| { let a = nx(&mut it).parse().unwrap(); let b = nx(&mut it).parse().unwrap(); let c = nx(&mut it).parse().unwrap(); [a, b, c] }).collect();
    assert_eq!(nx(&mut it), "G");
    let mut sign = [[1i8; NMAX]; M1];
    for i in 0..M1 { for t in 0..n { sign[i][t] = nx(&mut it).parse::<i8>().unwrap(); } }
    assert_eq!(nx(&mut it), "T"); let th = [nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap()];
    assert_eq!(nx(&mut it), "W"); let w = [nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap()];
    assert_eq!(nx(&mut it), "PD"); let pd: usize = nx(&mut it).parse().unwrap();
    assert_eq!(nx(&mut it), "H"); let k: usize = nx(&mut it).parse().unwrap();
    let hx = (0..k).map(|_| [nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap()]).collect();
    assert_eq!(nx(&mut it), "Q"); let k: usize = nx(&mut it).parse().unwrap();
    let gate = (0..k).map(|_| [nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap(), nx(&mut it).parse().unwrap()]).collect();
    Tv { n, samples, sign, th, w, pd, hx, gate }
}

fn build<T: Num>(tv: &Tv) -> Filter<T> {
    let mut f = Filter::<T>::new(tv.n, tv.sign);
    for s in &tv.samples { f.push(T::c(s[0]), T::c(s[1]), T::c(s[2])); }
    f.derive(); f
}

// xorshift: ismetelheto "szenzor"-sorozat az idomereshez
struct Xs(u32);
impl Xs { fn next(&mut self) -> f32 { let mut x = self.0; x ^= x << 13; x ^= x >> 17; x ^= x << 5; self.0 = x; (x as f32 / u32::MAX as f32) - 0.5 } }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).cloned().unwrap_or("tv20.txt".to_string());
    let iters: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(20000);
    let tv = read_tv(&path);
    println!("== ELLENORZES a Python (f64) referenciaval | N = {}", tv.n);
    let f: Filter<f32> = build(&tv);
    println!("theta Rust {:?}  Python {:?}  max eltérés {:.2e}", f.th, tv.th, (f.th[0] - tv.th[0]).abs().max((f.th[1] - tv.th[1]).abs()));
    println!("W     Rust {:?}  Python {:?}  max eltérés {:.2e}", f.w, tv.w, (0..3).map(|i| (f.w[i] - tv.w[i]).abs()).fold(0.0f32, f32::max));
    println!("nem-pozitiv-definit ellipszoidok: Rust {}  Python {}", f.n_np, tv.pd);
    let mut mx = 0.0f32; let mut mr = 0.0f32;
    for h in &tv.hx { let r = f.hx(h[0], h[1]); mx = mx.max((r - h[2]).abs()); mr = mr.max((r - h[2]).abs() / h[2].abs().max(1e-6)); }
    println!("burok-tartoszam (16 irany): max abs eltérés {:.2e}, max relatív {:.2e}", mx, mr);
    let (mut dmax, mut dsum, mut cnt, mut ndiff) = (0.0f32, 0.0f32, 0, 0);
    for g in &tv.gate { let (u, st) = f.gate(g[0], g[1]); if g[2].is_nan() { continue; } let d = (u - g[2]).abs(); dmax = dmax.max(d); dsum += d; cnt += 1; if d > 0.03 { ndiff += 1; } let _ = st; }
    println!("kapu (u): {} eset, max eltérés {:.4}, átlag {:.4}, >0.03: {} (a Python rács 0.025 lépésközű)", cnt, dmax, dsum / cnt as f32, ndiff);

    // inkrementalis vs teljes ujraszamolas: csuszo ablak 1000 lepes utan
    let mut rng = Xs(12345); let mut g: Filter<f32> = Filter::new(tv.n, tv.sign);
    let mut hist: Vec<[f32; 3]> = Vec::new();
    for t in 0..1000 { let v = 0.8 * (t as f32 * 0.07).sin() + 0.2 * rng.next(); let u = 1.0 + 0.5 * (t as f32 * 0.13).cos() + 0.2 * rng.next(); let y = 0.97 * v + 0.49 * u + 0.05 * rng.next(); g.push(v, u, y); hist.push([v, u, y]); }
    let mut refer: Filter<f32> = Filter::new(tv.n, tv.sign);
    for t in (1000 - tv.n)..1000 { refer.push(hist[t][0], hist[t][1], hist[t][2]); }   // az ablak sorrendje elter (slot) -> az elojel-oszlop ugyanaz a slot-hoz tartozik
    // a slot-hozzarendelest egyeztetjuk: g slotja (t mod n); refer slotja (t - (1000-n)) -> ugyanaz ha (1000-n) % n == 0, egyebkent a tabla ciklikus eltolasa
    let mut maxrel = 0.0f32; let shift = (1000 - tv.n) % tv.n;
    let mut refer2: Filter<f32> = Filter::new(tv.n, tv.sign);
    // g-ben a t. minta a (t mod n) oszlopot kapta: ujraepitjuk azonos oszlopokkal
    for t in (1000 - tv.n)..1000 { let s = t % tv.n; refer2.slot = s; refer2.push(hist[t][0], hist[t][1], hist[t][2]); }
    let _ = (&refer, shift);
    for i in 0..M1 { for k in 0..3 { let a = g.ai[i][k]; let b = refer2.ai[i][k]; maxrel = maxrel.max((a - b).abs() / b.abs().max(1.0)); } for k in 0..2 { let a = g.bi[i][k]; let b = refer2.bi[i][k]; maxrel = maxrel.max((a - b).abs() / b.abs().max(1.0)); } }
    println!("inkrementális (1000 csúszó lépés) vs. teljes újraszámolás, max relatív eltérés a B_i, A_i összegekben: {:.2e}", maxrel);

    println!("\n== MEMORIA (f32): Filter szerkezet {} bajt (ebbol az elojel-tabla {} bajt i8-kent, bitekkent ~{} bajt)", std::mem::size_of::<Filter<f32>>(), M1 * NMAX, (M1 * NMAX + 7) / 8);

    // ---- IDOMERES (PC, f32, egy mag): valodi ablak, valodi kapu-esetek
    let mut rng = Xs(777); let _ = rng.next();
    let mut flt: Filter<f32> = build(&tv);
    let mut acc = 0.0f32;
    let t0 = Instant::now();
    for k in 0..iters { let sm = tv.samples[k % tv.n]; flt.push(sm[0], sm[1], sm[2]); acc += flt.th[0]; }
    let t_push = t0.elapsed().as_secs_f64() / iters as f64;
    let t0 = Instant::now(); for _ in 0..iters { flt.derive(); acc += flt.th[1]; } let t_derive = t0.elapsed().as_secs_f64() / iters as f64;
    let mut gt: Vec<(f64, u8, f32, f32)> = Vec::new();
    for g in &tv.gate {
        let reps = 2000usize; let t0 = Instant::now(); let mut st = 0u8;
        for _ in 0..reps { let (u, s_) = flt.gate(g[0], g[1]); acc += u; st = s_; }
        gt.push((t0.elapsed().as_secs_f64() / reps as f64, st, g[0], g[1]));
    }
    let mut ts: Vec<f64> = gt.iter().map(|x| x.0).collect(); ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let st_hist = { let mut h = [0usize; 4]; for x in &gt { h[x.1 as usize] += 1; } h };
    println!("\n== IDOMERES (PC, WebAssembly/V8, f32, 1 mag; {} iteráció; ellenőrző összeg {:.3})", iters, acc);
    println!("push (inkrementális ablak-frissítés): {:.2} µs | derive (θ, W, 99 ellipszoid): {:.2} µs", t_push * 1e6, t_derive * 1e6);
    println!("kapu a {} valódi esetre: min {:.2} µs, medián {:.2} µs, max {:.2} µs | állapotok [átengedett, korlátozott, vész, tartalék] = {:?}", ts.len(), ts[0] * 1e6, ts[ts.len() / 2] * 1e6, ts[ts.len() - 1] * 1e6, st_hist);
    println!("egy vezérlési lépés (push + derive + kapu): medián {:.2} µs, legrosszabb {:.2} µs", (t_push + t_derive + ts[ts.len() / 2]) * 1e6, (t_push + t_derive + ts[ts.len() - 1]) * 1e6);

    // ---- MUVELETSZAMLALAS (ugyanaz a kod, Cn tipussal, ugyanazon a valodi ablakon)
    let mut c: Filter<Cn> = Filter::new(tv.n, tv.sign);
    for s_ in &tv.samples { c.push(Cn(s_[0]), Cn(s_[1]), Cn(s_[2])); }
    let snap = |a: [u64; 5], b: [u64; 5]| -> [u64; 5] { [b[0] - a[0], b[1] - a[1], b[2] - a[2], b[3] - a[3], b[4] - a[4]] };
    let sm = tv.samples[0];
    let a0 = counts(); c.push(Cn(sm[0]), Cn(sm[1]), Cn(sm[2])); let a1 = counts(); c.derive(); let a2 = counts();
    let dp = snap(a0, a1); let dd = snap(a1, a2);
    println!("\n== MUVELETSZAM / LEPES (add-sub | mul | div | sqrt | cmp)");
    println!("{:26} add {:6}  mul {:6}  div {:5}  sqrt {:5}  cmp {:6} | flop {:6}", "push", dp[0], dp[1], dp[2], dp[3], dp[4], dp[0] + dp[1]);
    println!("{:26} add {:6}  mul {:6}  div {:5}  sqrt {:5}  cmp {:6} | flop {:6}", "derive", dd[0], dd[1], dd[2], dd[3], dd[4], dd[0] + dd[1]);
    let mut gs: Vec<([u64; 5], u8)> = Vec::new();
    for g in &tv.gate { let b = counts(); let (_, st) = c.gate(Cn(g[0]), Cn(g[1])); gs.push((snap(b, counts()), st)); }
    gs.sort_by_key(|x| x.0[0] + x.0[1]);
    let show = |name: &str, x: &([u64; 5], u8)| println!("{:26} add {:6}  mul {:6}  div {:5}  sqrt {:5}  cmp {:6} | flop {:6}  (állapot {})", name, x.0[0], x.0[1], x.0[2], x.0[3], x.0[4], x.0[0] + x.0[1], x.1);
    show("kapu, legkisebb", &gs[0]); show("kapu, medián", &gs[gs.len() / 2]); show("kapu, LEGROSSZABB", &gs[gs.len() - 1]);
    let w = &gs[gs.len() - 1].0;
    println!("\nWCET-jellegű összeg (push + derive + a legrosszabb kapu): flop {} | div {} | sqrt {} | cmp {}", dp[0] + dp[1] + dd[0] + dd[1] + w[0] + w[1], dp[2] + dd[2] + w[2], dp[3] + dd[3] + w[3], dp[4] + dd[4] + w[4]);
}
