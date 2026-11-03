"""
E25: DRON, a Kolumban-halmaz a DEKLARALT TERVEZESI TARTOMANNYAL (elozetes) is, + a csak-tartomany (alkotmanyos) valtozat, + sertett tartomany. SZIMULACIO.
Ugyanaz a feladat, mint E24 (kolumban_dron.py). A deklarald tartomany (a drón specifikacioja): a in [0.95, 0.99], b in B_NOM*[0.7, 1.3]  (teher +-30%).
Uj modszerek:
  DOBOZ  : CSAK a deklarald tartomany (nincs adat, nincs becsles): a szuro a doboz 4 csucsara a legrosszabb esetet nezi  -> determinisztikus, O(1); biztonsagos, HA a tartomany igaz
  RACS+D / IRANY+D / BUROK+D : a Kolumban-halmaz METSZVE a dobozzal (a burok-halmazt 64 irany tartoszamaval kulso sokszogge lapitjuk, majd a dobozzal vagjuk -> tartalmazza a halmaz es a doboz metszetet)
         ures metszet = az adat ellentmond a deklarald tartomanynak -> tartalek: a doboz (es jelezzuk)
Ket kornyezet:  ervenyes : a valodi b a deklarald +-30%-on belul;   serul : a valodi b a +-45%-ban szorodik (a deklarald +-30% HAMIS az ablakok ~1/3-aban).
Mertek (E24 + ures metszet aranya; a "serul" kornyezetben kulon az ablakokra, ahol a valodi theta a dobozon KIVUL van).
"""
import sys, time
import numpy as np
from scipy.stats import chi2, binomtest
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from kolumban_arx import build, ellipsoids, support_q, q, m
from kolumban_vezerles_futasido import grid_members, radial, wilson
import kolumban_dron as D

DT, G, GDT, B_NOM, A_NOM, V, UG, K = D.DT, D.G, D.GDT, D.B_NOM, D.A_NOM, D.V, D.UG, D.K
C95 = D.C95
A_LO, A_HI = 0.95, 0.99; B_LO, B_HI = 0.7 * B_NOM, 1.3 * B_NOM
BOX = np.array([[A_LO, B_LO], [A_HI, B_LO], [A_HI, B_HI], [A_LO, B_HI]])
NAMES = ["NONE", "NOM", "CE", "ASY", "EKF-nom", "EKF-becs", "DOBOZ", "RACS", "IRANY", "BUROK", "RACS+D", "IRANY+D", "BUROK+D", "ORA"]
ANG64 = np.linspace(0, 2 * np.pi, 64, endpoint=False); C64 = np.stack([np.cos(ANG64), np.sin(ANG64)], 1)
HALF_BOX = [(np.array([1.0, 0.0]), A_HI), (np.array([-1.0, 0.0]), -A_LO), (np.array([0.0, 1.0]), B_HI), (np.array([0.0, -1.0]), -B_LO)]


def clip(poly, n, d):
    out = []; L = len(poly)
    for i in range(L):
        p, q_ = poly[i], poly[(i + 1) % L]; fp, fq = n @ p - d, n @ q_ - d
        if fp <= 0: out.append(p)
        if fp * fq < 0: out.append(p + fp / (fp - fq) * (q_ - p))
    return np.array(out)


def window(N, noise, env, rng):
    dl = 0.3 if env == "ervenyes" else 0.45
    a0 = rng.uniform(A_LO, A_HI); b0 = B_NOM * (1 + rng.uniform(-dl, dl)); th0 = np.array([a0, b0])
    if noise == "laplace": w = rng.laplace(size=N) * 0.05 / np.sqrt(2)
    elif noise == "szeles": w = rng.laplace(size=N) * 0.12 / np.sqrt(2)
    else: w = rng.standard_t(1.5, size=N) * 0.03
    d = rng.normal(size=N) * D.SD_D
    v = np.zeros(N + 1); u = np.zeros(N); v[0] = rng.normal() * 0.3
    for t in range(N):
        u[t] = np.clip(1.0 - D.KFB * v[t] + d[t], 0.0, 3.0); v[t + 1] = a0 * v[t] + b0 * u[t] - GDT + w[t]
    return np.stack([v[:-1], u], 1), v[1:] + GDT, th0


def run(N, noise, env, R, rng):
    thn = np.array([A_NOM, B_NOM]); ph = np.linspace(0, 2 * np.pi, 16, endpoint=False); dirs = np.stack([np.cos(ph), np.sin(ph)], 1)
    S = dict(vany=np.zeros((R, len(NAMES)), bool), viol=np.zeros((R, len(NAMES))), inf=np.zeros((R, len(NAMES))), dev=[[] for _ in NAMES],
             cov={k: 0 for k in ("EKF-becs", "RACS", "IRANY", "BUROK", "RACS+D", "IRANY+D", "BUROK+D")}, empty={k: 0 for k in ("RACS+D", "IRANY+D", "BUROK+D")},
             empty_out={k: 0 for k in ("RACS+D", "IRANY+D", "BUROK+D")}, n_out=0, outbox=np.zeros(R, bool), t_set={k: 0.0 for k in ("RACS", "IRANY", "BUROK", "BUROK+D", "IRANY+D")},
             t_q={nm: 0.0 for nm in NAMES}, nq=0, skipped=0)
    i = 0
    while i < R:
        Phi, y, th0 = window(N, noise, env, rng)
        sps = build(Phi, y, "sign", rng); th = sps.theta_ls; W = sps.W; E = ellipsoids(sps)
        if (~E['pd']).sum() >= q: S['skipped'] += 1; continue
        s2 = ((y - Phi @ th) ** 2).sum() / (N - 2); Ri = np.linalg.inv(Phi.T @ Phi); se = np.sqrt(np.diag(Ri) * s2)
        t = time.perf_counter(); Mg, ntg = grid_members(sps, th, se); tg = time.perf_counter() - t
        if len(Mg) == 0: S['skipped'] += 1; continue
        t = time.perf_counter(); pts, nti = radial(sps, th, dirs, W); ti = time.perf_counter() - t
        t = time.perf_counter(); E = ellipsoids(sps); tu = time.perf_counter() - t
        def hx(Ct):
            Cz = Ct @ W; nrm = np.linalg.norm(Cz, axis=1) + 1e-300
            return support_q(E, Cz / nrm[:, None]) * nrm
        # --- deklarald tartomannyal metszett halmazok (csucslistak) ---
        t = time.perf_counter(); s64 = C64 @ th + hx(C64); poly = BOX.copy()
        for c, s_ in zip(C64, s64):
            poly = clip(poly, c, s_)
            if len(poly) < 3: break
        e_b = len(poly) < 3; Vb = BOX if e_b else poly; tbd = tu + time.perf_counter() - t
        t = time.perf_counter()
        try: hp = pts[ConvexHull(pts).vertices]
        except Exception: hp = pts
        poly_i = hp
        for n_, d_ in HALF_BOX:
            poly_i = clip(poly_i, n_, d_)
            if len(poly_i) < 3: break
        e_i = len(poly_i) < 3; Vi = BOX if e_i else poly_i; tid = ti + time.perf_counter() - t
        inb = (Mg[:, 0] >= A_LO) & (Mg[:, 0] <= A_HI) & (Mg[:, 1] >= B_LO) & (Mg[:, 1] <= B_HI)
        e_r = not inb.any(); Vr = BOX if e_r else Mg[inb]
        for k_, v_ in (("RACS", tg), ("IRANY", ti), ("BUROK", tu), ("BUROK+D", tbd), ("IRANY+D", tid)): S['t_set'][k_] += v_
        outb = not ((A_LO <= th0[0] <= A_HI) and (B_LO <= th0[1] <= B_HI)); S['outbox'][i] = outb; S['n_out'] += int(outb)
        for k_, e_ in (("RACS+D", e_r), ("IRANY+D", e_i), ("BUROK+D", e_b)):
            S['empty'][k_] += int(e_)
            if outb: S['empty_out'][k_] += int(e_)
        thE1, PE1 = D.rls(Phi, y, D.R_NOM); thE2, PE2 = D.rls(Phi, y, s2)
        ell = lambda Pm, t_, Sg: np.abs(Pm @ t_ - GDT) + np.sqrt(C95 * np.einsum('ij,jl,il->i', Pm, Sg, Pm))
        vw = lambda Vv: (lambda Pm: np.abs(Pm @ Vv.T - GDT).max(1))
        fn = {"NONE": None, "NOM": lambda Pm: np.abs(Pm @ thn - GDT), "CE": lambda Pm: np.abs(Pm @ th - GDT),
              "ASY": lambda Pm: ell(Pm, th, s2 * Ri), "EKF-nom": lambda Pm: ell(Pm, thE1, PE1), "EKF-becs": lambda Pm: ell(Pm, thE2, PE2),
              "DOBOZ": vw(BOX), "RACS": lambda Pm: np.abs(Pm @ Mg.T - GDT).max(1), "IRANY": lambda Pm: np.abs(Pm @ pts.T - GDT).max(1),
              "BUROK": lambda Pm: np.maximum(Pm @ th + hx(Pm) - GDT, -(Pm @ th) + hx(-Pm) + GDT),
              "RACS+D": vw(Vr), "IRANY+D": vw(Vi), "BUROK+D": vw(Vb), "ORA": lambda Pm: np.abs(Pm @ th0 - GDT)}
        dd = th0 - thE2; S['cov']["EKF-becs"] += int(dd @ np.linalg.solve(PE2, dd) <= C95)
        S['cov']["RACS"] += int(sps.member(th0[None])[0]); S['cov']["IRANY"] += int(Path(pts).contains_point(th0))
        ang = np.linspace(0, 2 * np.pi, 720, endpoint=False); C720 = np.stack([np.cos(ang), np.sin(ang)], 1)
        z0 = C720 @ np.linalg.solve(W, th0 - th); S['cov']["BUROK"] += int((z0 <= support_q(E, C720) + 1e-9).all())
        S['cov']["RACS+D"] += int(S['cov'].get("_", 0) * 0 + (sps.member(th0[None])[0] and not outb))
        S['cov']["IRANY+D"] += int(Path(pts).contains_point(th0) and not outb)
        S['cov']["BUROK+D"] += int(bool((z0 <= support_q(E, C720) + 1e-9).all()) and not outb)
        xs = rng.uniform(-1, 1, K); un = rng.uniform(0.2, 1.8, K)
        for j, nm in enumerate(NAMES):
            if nm == "NONE": u_ = un.copy()
            else:
                t = time.perf_counter(); u_ = D.robust_u(fn[nm], xs, un); S['t_q'][nm] += time.perf_counter() - t
            nanm = np.isnan(u_); ok = ~nanm
            vv = (np.abs(th0[0] * xs[ok] + th0[1] * u_[ok] - GDT) > V + 1e-12) if ok.any() else np.zeros(0, bool)
            S['viol'][i, j] = vv.mean() if ok.any() else 0.0; S['vany'][i, j] = vv.any() if ok.any() else False
            S['inf'][i, j] = nanm.mean(); S['dev'][j] += list(np.abs(u_[ok] - un[ok]))
        S['nq'] += K; i += 1
    return S


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    envs = sys.argv[2].split(",") if len(sys.argv) > 2 else ["ervenyes", "serul"]
    conds = sys.argv[3].split(",") if len(sys.argv) > 3 else ["laplace", "szeles", "heavy"]
    Ns = [int(v) for v in sys.argv[4].split(",")] if len(sys.argv) > 4 else [12, 20, 40]
    print(f"SZIMULACIO (E25). deklarald tartomany: a in [{A_LO},{A_HI}], b in [{B_LO:.3f},{B_HI:.3f}] (nominalis +-30%); {R} ablak/cella, V={V} m/s")
    store = {}
    for env in envs:
        for noise in conds:
            for N in Ns:
                S = run(N, noise, env, R, np.random.default_rng(40000 + N + 11 * len(noise) + 7 * len(env))); store[(env, noise, N)] = S
                nout = S['n_out']
                print(f"\n=== kornyezet: {env}, zaj: {noise}, N={N} | ablak {R} (vakuum kihagyva {S['skipped']}) | valodi theta a dobozon KIVUL: {nout} ablak | ms/ablak halmaz: racs {S['t_set']['RACS']/R*1e3:.1f}, irany {S['t_set']['IRANY']/R*1e3:.1f}, burok {S['t_set']['BUROK']/R*1e3:.2f}, burok+D {S['t_set']['BUROK+D']/R*1e3:.2f}, irany+D {S['t_set']['IRANY+D']/R*1e3:.1f}")
                print(f"{'modszer':10s}{'serulo dontes %':>16s}{'van serules ablak % [Wilson]':>32s}{'nincs megold. %':>17s}{'atl |u-unom|':>14s}{'lek. us/dontes':>16s}{'lefedes':>9s}{'ures metszet %':>16s}{'(dobozon kivul: ures %)':>25s}")
                for j, nm in enumerate(NAMES):
                    va = S['vany'][:, j].sum(); lo, hi = wilson(int(va), R)
                    cv = f"{S['cov'][nm]/R:9.3f}" if nm in S['cov'] else f"{'-':>9s}"
                    em = f"{S['empty'][nm]/R*100:16.1f}" if nm in S['empty'] else f"{'-':>16s}"
                    emo = f"{S['empty_out'][nm]/max(nout,1)*100:25.1f}" if (nm in S['empty_out'] and nout) else f"{'-':>25s}"
                    print(f"{nm:10s}{S['viol'][:, j].mean()*100:16.3f}{va/R*100:14.1f} [{lo*100:.1f},{hi*100:.1f}]{S['inf'][:, j].mean()*100:17.1f}{np.mean(S['dev'][j]) if S['dev'][j] else float('nan'):14.3f}{(S['t_q'][nm]/S['nq']*1e6 if nm != 'NONE' else 0):16.0f}{cv}{em}{emo}", flush=True)
    print("\n=== PAROSITOTT TESZT (ablak-szintu 'van serules', McNemar pontos binomialis) ===")
    pairs = [("BUROK+D", "EKF-becs"), ("BUROK+D", "BUROK"), ("BUROK+D", "DOBOZ"), ("BUROK+D", "RACS+D"), ("IRANY+D", "RACS+D"), ("EKF-becs", "DOBOZ")]
    for (env, noise, N), S in store.items():
        print(f"\n{env}, {noise}, N={N}:")
        for a, b in pairs:
            ja, jb = NAMES.index(a), NAMES.index(b); x, y_ = S['vany'][:, ja], S['vany'][:, jb]
            b01 = int((x & ~y_).sum()); b10 = int((~x & y_).sum()); p = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 > 0 else 1.0
            print(f"  {a:8s} vs {b:9s}: csak {a} serul {b01:4d}, csak {b} serul {b10:4d} | p = {p:.4f}" + ("  *" if p < 0.05 else ""))
