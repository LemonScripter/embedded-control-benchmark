"""
E26: ZART HURKU REPULES (drón függőleges tengely), architektúrák összevetése. SZIMULACIO.
Epizód: T=400 lépés (20 s, dt=0.05). A pilóta/autopilot sebesség-referenciát követ (PI + előrecsatolás a NOMINALIS modellel), a referencia szakaszonként állandó, +-1.4 m/s (a biztonsági határ |v+|<=1 m/s, ezért a pilóta olykor többet kér, mint amit szabad).
A 120. lépésben TEHERVALTAS: a tolóerő-hatásfok b a névleges (1+d)-szorosa lesz (d ~ U(-0.3,0.3), érvényes tartomány; vagy U(-0.45,0.45), sérülő tartomány), a légellenállás a is új értéket kap.
Az első Nw lépés felszállás-előtti azonosítás (nincs szűrő, nem számít a mérőszámba). Utána minden architektúra a saját adatából (utolsó Nw minta, csúszó ablak) ítéli meg a drónt, és a szűrő a pilóta u_nom parancsát
a legközelebbi biztonságosra módosítja: |v+|<=V MINDEN hihető (a,b)-re. Ha nincs biztonságos parancs: vész-parancs = a legkisebb legrosszabb-esetű (és számoljuk).
Architektúrák (ugyanaz a szél-, dither-, referencia-sorozat mindegyiknek, párosítva):
  PID       : nincs szűrő (a pilóta parancsa megy a motorra)
  NOM       : a szűrő fix névleges modellel (nincs tanulás)
  EKF-nom   : parameter-EKF (RLS + parameter-random-walk), tervezői zajbecslés R=0.05^2, helyes előzetes (névleges +- szórás), 95% ellipszoid
  EKF-becs  : ugyanaz, R a maradékból (csúszó)
  BUROK     : Kolumbán-burok (csúszó Nw ablak, m=100, q=5), zajmodell és előzetes NELKÜL
  BUROK+D   : burok metszve a DEKLARALT tartománnyal (a, b +-30%), ures metszet -> tartalék: csak az adatalapu burok (E25b)
  DOBOZ     : csak a deklarald tartomány (nincs adat)
  ORA       : a valodi (a,b) ismert (nem megvalosithato; a felso korlat)
Mertek: serulo lepes % (|v+|>V a valodi dronra), "van serules" epizod %, beavatkozas |u-u_nom|, kovetesi hiba RMSE(v - clip(v_ref,+-V)), vesz-lepesek %, ms/lepes (PC).
"""
import sys, time
import numpy as np
from scipy.stats import chi2, binomtest
from kolumban_arx import build, ellipsoids, support_q, q, m
import kolumban_dron as D
import kolumban_dron_elozetes as E

DT, G, GDT, B_NOM, A_NOM, V = D.DT, D.G, D.GDT, D.B_NOM, D.A_NOM, D.V
UG = D.UG; C95 = D.C95; BOX = E.BOX
T, TCH = 400, 120
KP, KI = 0.8, 0.4
SD_D = 0.2
ARCH = ["PID", "NOM", "EKF-nom", "EKF-becs", "BUROK", "BUROK+D", "DOBOZ", "ORA"]
P0 = np.diag([0.02 ** 2, 0.10 ** 2]); QTH = np.diag([1e-7, 1e-4]); R_NOM = 0.05 ** 2
C64 = E.C64


def noise_seq(noise, rng):
    if noise == "laplace": return rng.laplace(size=T) * 0.05 / np.sqrt(2)
    if noise == "szeles": return rng.laplace(size=T) * 0.12 / np.sqrt(2)
    return rng.standard_t(1.5, size=T) * 0.03


def make_episode(noise, env, rng):
    dl = 0.3 if env == "ervenyes" else 0.45
    a1 = rng.uniform(0.95, 0.99); a2 = rng.uniform(0.95, 0.99); b1 = B_NOM * (1 + rng.uniform(-0.1, 0.1)); b2 = B_NOM * (1 + rng.uniform(-dl, dl))
    th = np.array([[a1, b1]] * TCH + [[a2, b2]] * (T - TCH))
    vref = np.zeros(T); t = 20
    while t < T:
        L = int(rng.integers(20, 50)); vref[t:t + L] = rng.uniform(-1.4, 1.4); t += L
    return dict(th=th, w=noise_seq(noise, rng), d=rng.normal(size=T) * SD_D, vref=vref)


def pid_u(v, vr, integ):
    e = vr - v
    return (vr - A_NOM * v + GDT) / B_NOM + KP * e + KI * integ


class Est:
    """parameter-EKF (RLS + parameter-random-walk)"""
    def __init__(self, adapt_R, robust=False):
        self.th = np.array([A_NOM, B_NOM]); self.P = P0.copy(); self.R = R_NOM; self.adapt = adapt_R; self.rob = robust
    def update(self, phi, y):
        self.P = self.P + QTH
        inn = y - phi @ self.th; S = phi @ self.P @ phi + self.R
        if self.rob: inn = float(np.clip(inn, -2.5 * np.sqrt(S), 2.5 * np.sqrt(S)))   # Huber-szeru vagas (robusztus EKF)
        K = self.P @ phi / S
        self.th = self.th + K * inn; self.P = self.P - np.outer(K, phi @ self.P)
        if self.adapt: self.R = max(1e-4, 0.97 * self.R + 0.03 * inn ** 2)


def run_arch(name, ep, Nw, rng_k, tm, trace=None):
    th_t, w, d, vref = ep["th"], ep["w"], ep["d"], ep["vref"]
    v = 0.0; integ = 0.0; Phi_h, y_h = [], []
    est = Est(False) if name == "EKF-nom" else (Est(True) if name == "EKF-becs" else (Est(True, True) if name == "EKF-rob" else None))
    viol = np.zeros(T, bool); dev = np.zeros(T); trk = np.zeros(T); emerg = np.zeros(T, bool); vs = np.zeros(T); noexc = np.zeros(T, bool)
    for t in range(T):
        e = vref[t] - v; integ = float(np.clip(integ + e * DT, -2, 2))
        u_nom = float(np.clip(pid_u(v, vref[t], integ) + d[t], 0.0, 3.0))
        active = (t >= Nw) and name != "PID"
        u = u_nom; bi = (float('nan'), float('nan'))
        if active:
            Pm = np.stack([np.full_like(UG, v), UG], 1)
            if name == "NOM": worst = np.abs(Pm @ np.array([A_NOM, B_NOM]) - GDT)
            elif name == "ORA": worst = np.abs(Pm @ th_t[t] - GDT)
            elif name == "DOBOZ": worst = np.abs(Pm @ BOX.T - GDT).max(1); bi = (float(BOX[:, 1].min()), float(BOX[:, 1].max()))
            elif name in ("EKF-nom", "EKF-becs", "EKF-rob"):
                t0 = time.perf_counter(); thh, PP = est.th, est.P
                worst = np.abs(Pm @ thh - GDT) + np.sqrt(C95 * np.einsum('ij,jl,il->i', Pm, PP, Pm)); tm[name] += time.perf_counter() - t0
                bi = (float(thh[1] - np.sqrt(C95 * PP[1, 1])), float(thh[1] + np.sqrt(C95 * PP[1, 1])))
            else:   # BUROK / BUROK+D
                t0 = time.perf_counter()
                Phi = np.array(Phi_h[-Nw:]); yy = np.array(y_h[-Nw:])
                Rm = Phi.T @ Phi / Nw
                if np.linalg.eigvalsh(Rm)[0] < 1e-7:        # nincs elegendo gerjesztes: az ablak szingularis, adatalapu halmaz nem kepezheto
                    vac = True; noexc[t] = True
                else:
                    sps = build(Phi, yy, "sign", rng_k); thl = sps.theta_ls; W = sps.W; Eel = ellipsoids(sps)
                    vac = (~Eel['pd']).sum() >= q
                if vac:
                    worst = np.abs(Pm @ np.array([A_NOM, B_NOM]) - GDT) if name == "BUROK" else np.abs(Pm @ BOX.T - GDT).max(1)
                else:
                    def hx(Ct):
                        Cz = Ct @ W; nrm = np.linalg.norm(Cz, axis=1) + 1e-300
                        return support_q(Eel, Cz / nrm[:, None]) * nrm
                    w_hull = np.maximum(Pm @ thl + hx(Pm) - GDT, -(Pm @ thl) + hx(-Pm) + GDT)
                    bi = (float(thl[1] - hx(np.array([[0.0, -1.0]]))[0]), float(thl[1] + hx(np.array([[0.0, 1.0]]))[0]))
                    if name == "BUROK": worst = w_hull
                    else:
                        poly = BOX.copy(); s64 = C64 @ thl + hx(C64)
                        for c, s_ in zip(C64, s64):
                            poly = E.clip(poly, c, s_)
                            if len(poly) < 3: break
                        worst = w_hull if len(poly) < 3 else np.abs(Pm @ poly.T - GDT).max(1)
                        if len(poly) >= 3: bi = (float(poly[:, 1].min()), float(poly[:, 1].max()))
                tm[name] += time.perf_counter() - t0
            ok = worst <= V
            if ok.any(): u = float(UG[np.argmin(np.where(ok, np.abs(UG - u_nom), np.inf))])
            else: u = float(UG[np.argmin(worst)]); emerg[t] = True
        a_t, b_t = th_t[t]
        v_next = a_t * v + b_t * u - GDT + w[t]
        viol[t] = abs(a_t * v + b_t * u - GDT) > V + 1e-12; dev[t] = abs(u - u_nom); trk[t] = (v_next - np.clip(vref[t], -V, V)) ** 2; vs[t] = abs(a_t * v + b_t * u - GDT)
        if trace is not None: trace.append((t, float(v), float(u_nom), float(u), bool(viol[t]), bool(emerg[t]), bi[0], bi[1], float(b_t)))
        phi = np.array([v, u]); y = v_next + GDT
        if est is not None: t0 = time.perf_counter(); est.update(phi, y); tm[name] += time.perf_counter() - t0
        Phi_h.append(phi); y_h.append(y); v = v_next
    sl = slice(Nw, T); post = slice(TCH, TCH + Nw + 20)
    return dict(viol=viol[sl].mean(), anyv=bool(viol[sl].any()), post=viol[post].mean(), dev=dev[sl].mean(), trk=float(np.sqrt(trk[sl].mean())), emerg=emerg[sl].mean(), vmax=float(vs[sl].max()), noexc=noexc[sl].mean())


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    Nw = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    env = sys.argv[3] if len(sys.argv) > 3 else "ervenyes"
    noises = sys.argv[4].split(",") if len(sys.argv) > 4 else ["laplace", "szeles"]
    print(f"SZIMULACIO (E26, zart hurku repules). {R} epizod/cella, T={T} lepes, teherváltás a {TCH}. lepesben, csuszó ablak Nw={Nw}, kornyezet: {env}")
    for noise in noises:
        rng = np.random.default_rng(60000 + Nw + 7 * len(noise) + 3 * len(env)); res = {a: [] for a in ARCH}; tm = {a: 0.0 for a in ARCH}
        for _ in range(R):
            ep = make_episode(noise, env, rng)
            for a in ARCH: res[a].append(run_arch(a, ep, Nw, np.random.default_rng(int(rng.integers(1 << 30))), tm))
        print(f"\n=== zaj: {noise}, Nw={Nw}, kornyezet: {env} | epizod {R}")
        print(f"{'architektura':10s}{'serulo lepes %':>15s}{'van serules epizod % [Wilson]':>34s}{'teherv. utan %':>16s}{'beavatk.':>10s}{'kov. hiba':>10s}{'vesz %':>8s}{'nincs ger.%':>12s}{'ms/lepes':>10s}")
        for a in ARCH:
            r = res[a]; va = sum(x['anyv'] for x in r); lo, hi = wilson(va, R)
            print(f"{a:10s}{np.mean([x['viol'] for x in r])*100:15.3f}{va/R*100:14.1f} [{lo*100:.1f},{hi*100:.1f}]{np.mean([x['post'] for x in r])*100:16.2f}{np.mean([x['dev'] for x in r]):10.3f}{np.mean([x['trk'] for x in r]):10.3f}{np.mean([x['emerg'] for x in r])*100:8.2f}{np.mean([x['noexc'] for x in r])*100:12.2f}{tm[a]/(R*(T-Nw))*1e3:10.3f}", flush=True)
        print("  PAROSITOTT (epizod-szintu 'van serules', McNemar):")
        for x_, y_ in (("BUROK", "EKF-becs"), ("BUROK", "EKF-nom"), ("BUROK", "NOM"), ("BUROK+D", "EKF-becs"), ("BUROK+D", "BUROK"), ("BUROK+D", "DOBOZ"), ("EKF-becs", "PID")):
            a1 = np.array([z['anyv'] for z in res[x_]]); a2 = np.array([z['anyv'] for z in res[y_]]); b01 = int((a1 & ~a2).sum()); b10 = int((~a1 & a2).sum())
            p = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 else 1.0
            print(f"    {x_:8s} vs {y_:9s}: csak {x_} serul {b01:4d}, csak {y_} serul {b10:4d} | p = {p:.4f}" + ("  *" if p < 0.05 else ""))
