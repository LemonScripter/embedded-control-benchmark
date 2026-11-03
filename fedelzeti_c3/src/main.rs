//! ESP32-C3 (RV32IMC, FPU nelkul) cel: ugyanaz a Kolumban-burok kapu algoritmus, #![no_std], f32 SZOFTVERES lebegopontos muveletekkel.
//! A gazdagep (Python + Unicorn emulator) a memoriaba irja az elojel-tablat (SIGN), meghivja a c3_* fuggvenyeket, es utasitasokat szamol.
#![no_std]
#![no_main]
use core::ops::{Add, Div, Mul, Neg, Sub};
use core::panic::PanicInfo;

#[panic_handler]
fn panic(_: &PanicInfo) -> ! { loop {} }

core::arch::global_asm!(".section .text.init", ".global _start", "_start:", "la sp, _stack_top", "1: j 1b");

const M1: usize = 99; // m - 1
const NMAX: usize = 40;
const Q: usize = 5;
const GDT: f32 = 0.4905;
const VLIM: f32 = 1.0;
const A_NOM: f32 = 0.97;
const B_NOM: f32 = 0.4905;
const PD_EPS: f32 = 1e-6;

/// libm-mentes sqrt: bit-hack kezdoertek + 3 Newton-lepes (hiba < 1e-7 relativ); x <= 0 -> 0
fn sqrtf_soft(x: f32) -> f32 {
    if !(x > 0.0) { return 0.0; }
    if x == f32::INFINITY { return x; }
    let mut y = f32::from_bits((x.to_bits() >> 1) + 0x1FBD_1DF5);
    y = 0.5 * (y + x / y); y = 0.5 * (y + x / y); y = 0.5 * (y + x / y);
    y
}

trait Num: Copy + PartialOrd + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self> {
    fn c(x: f32) -> Self;
    fn sqrt(self) -> Self;
}
impl Num for f32 {
    fn c(x: f32) -> f32 { x }
    fn sqrt(self) -> f32 { sqrtf_soft(self) }
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

#[no_mangle] pub static mut SIGN: [[i8; NMAX]; M1] = [[0; NMAX]; M1];
#[no_mangle] pub static mut OUT: [f32; 8] = [0.0; 8];
#[no_mangle] pub static mut STATE: u32 = 0;
static mut FLT: Filter<f32> = Filter {
    n: 0, sign: [[0; NMAX]; M1], sv: [0.0; NMAX], su: [0.0; NMAX], sy: [0.0; NMAX], slot: 0, r: [0.0; 3], b0: [0.0; 2], bi: [[0.0; 2]; M1], ai: [[0.0; 3]; M1],
    th: [0.0; 2], w: [0.0; 3], zeta: [[0.0; 2]; M1], pinv: [[0.0; 3]; M1], kappa: [0.0; M1], pd: [false; M1], n_np: M1, ok: false };

#[no_mangle] pub extern "C" fn c3_init(n: u32) { unsafe { let f = &mut *core::ptr::addr_of_mut!(FLT); *f = Filter::new(n as usize, *core::ptr::addr_of!(SIGN)); } }
#[no_mangle] pub extern "C" fn c3_push(v: f32, u: f32, y: f32) { unsafe { (*core::ptr::addr_of_mut!(FLT)).push(v, u, y); } }
#[no_mangle] pub extern "C" fn c3_derive() { unsafe { (*core::ptr::addr_of_mut!(FLT)).derive(); } }
#[no_mangle] pub extern "C" fn c3_gate(v: f32, u_nom: f32) -> f32 { unsafe { let (u, s) = (*core::ptr::addr_of!(FLT)).gate(v, u_nom); STATE = s as u32; u } }
#[no_mangle] pub extern "C" fn c3_hx(pv: f32, pu: f32) -> f32 { unsafe { (*core::ptr::addr_of!(FLT)).hx(pv, pu) } }
#[no_mangle] pub extern "C" fn c3_dump() { unsafe { let f = &*core::ptr::addr_of!(FLT); OUT = [f.th[0], f.th[1], f.w[0], f.w[1], f.w[2], f.n_np as f32, 0.0, 0.0]; } }
