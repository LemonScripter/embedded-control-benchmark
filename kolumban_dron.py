"""
E24: DRON, fuggoleges tengely (landolas / tehermentes-teherrel valtas): biztonsagos repules ismeretlen dinamika es keves adat mellett. SZIMULACIO.
Fizika (20 Hz, dt=0.05 s): v+ = a v + b u - g dt + w ; v = fuggoleges sebesseg (m/s, felfele pozitiv), u = normalizalt tolóero-parancs (1 = hover a NOMINALIS teherrel),
  a = 1 - suruseg*dt (ismeretlen legellenallas, 0.95..0.99), b = tolóero-hatasfok (ISMERETLEN: a teher / akku miatt a nominalis +-30%), w = szel/turbulencia (szimmetrikus, fuggetlen; Laplace vagy nehez farku).
Adat: a teher valtasa utan csak N=12..40 minta (0.6..2 s) all rendelkezesre, kozben a hover-szabalyozo (u = 1 - k v + dither) repul (zart hurok, gyenge gerjesztes) -> INFORMACIOSZEGENY.
Biztonsagi szabaly (landolasi kaput / lombkoz-limit): a kovetkezo lepesben |v+| <= V MINDEN hihetonek tartott (a, b)-re. A pilota u_nom-ot keri; a szuro a legkozelebbi biztonsagos u-t adja, vagy "nincs biztonsagos u" (vegso eset: vesz-szabaly).
Biztonsagi hiba: a VALODI (a0, b0)-ra |v+| > V.
Modszerek (ugyanabbol az ablakbol, parositva):
  NONE : nincs szuro                      NOM  : fix nominalis modell (nincs tanulas)       CE : csak a legkisebb negyzetes theta^
  ASY  : klasszikus aszimptotikus 95% ellipszoid (normalitas)
  EKF-nom / EKF-becs : rekurziv (Kalman) parameter-becslo HELYES elozetessel (a nominalis +-szoras), R = tervezoi / a maradekbol; 95% chi2 ellipszoid a P-bol
  RACS : pontos SPS-halmaz racson (121x121)   IRANY : 16 irany   BUROK : lapitas (0 teszt)      [mind elozetes NELKUL, zajmodell NELKUL]
  ORA  : a valodi (a0, b0)
Mertek: serulo dontes %, ablak-szintu "van serules" (McNemar-parositott teszt), nincs megoldas %, atl |u-unom| (a kenyszer ara), a valodi theta lefedese, ido.
"""
import sys, time
import numpy as np
from scipy.stats import chi2, binomtest
from matplotlib.path import Path
from kolumban_arx import build, ellipsoids, support_q, q, m
from kolumban_vezerles_futasido import grid_members, radial, wilson

DT = 0.05; G = 9.81; GDT = G * DT; B_NOM = GDT; A_NOM = 0.97
V = 1.0
UG = np.linspace(0.0, 3.0, 121)
K = 40; KFB = 1.0; SD_D = 0.25
C95 = chi2.ppf(0.95, 2)
P0 = np.diag([0.02 ** 2, 0.10 ** 2]); R_NOM = 0.05 ** 2
NAMES = ["NONE", "NOM", "CE", "ASY", "EKF-nom", "EKF-becs", "RACS", "IRANY", "BUROK", "ORA"]


def window(N, noise, rng):
    a0 = rng.uniform(0.95, 0.99); b0 = B_NOM * (1 + rng.uniform(-0.3, 0.3)); th0 = np.array([a0, b0])
    if noise == "laplace": w = rng.laplace(size=N) * 0.05 / np.sqrt(2)
    elif noise == "szeles": w = rng.laplace(size=N) * 0.12 / np.sqrt(2)          # szeles nap: a tervezo R-je (0.05) alulmeretezett
    else: w = rng.standard_t(1.5, size=N) * 0.03
    d = rng.normal(size=N) * SD_D
    v = np.zeros(N + 1); u = np.zeros(N); v[0] = rng.normal() * 0.3
    for t in range(N):
        u[t] = np.clip(1.0 - KFB * v[t] + d[t], 0.0, 3.0); v[t + 1] = a0 * v[t] + b0 * u[t] - GDT + w[t]
    return np.stack([v[:-1], u], 1), v[1:] + GDT, th0


def robust_u(worst, xs, unom):
    out = np.full(len(xs), np.nan)
    for i, (x, un) in enumerate(zip(xs, unom)):
        Phi = np.stack([np.full_like(UG, x), UG], 1); ok = worst(Phi) <= V
        if ok.any(): out[i] = UG[np.argmin(np.where(ok, np.abs(UG - un), np.inf))]
    return out


def rls(Phi, y, R):
    th = np.array([A_NOM, B_NOM]); P = P0.copy()
    for f, yy in zip(Phi, y):
        k = P @ f / (f @ P @ f + R); th = th + k * (yy - f @ th); P = P - np.outer(k, f @ P)
    return th, P


def run(N, noise, R, rng):
    ang = np.linspace(0, 2 * np.pi, 720, endpoint=False); C720 = np.stack([np.cos(ang), np.sin(ang)], 1)
    ph = np.linspace(0, 2 * np.pi, 16, endpoint=False); dirs = np.stack([np.cos(ph), np.sin(ph)], 1)
    thn = np.array([A_NOM, B_NOM])
    S = dict(n=0, viol=np.zeros((R, len(NAMES))), vany=np.zeros((R, len(NAMES)), bool), inf=np.zeros((R, len(NAMES))), dev=[[] for _ in NAMES],
             cov={k: 0 for k in ("ASY", "EKF-nom", "EKF-becs", "RACS", "IRANY", "BUROK")}, t_set={"RACS": 0.0, "IRANY": 0.0, "BUROK": 0.0}, nt={"RACS": 0, "IRANY": 0},
             t_q={nm: 0.0 for nm in NAMES}, nq=0, skipped=0)
    i = 0
    while i < R:
        Phi, y, th0 = window(N, noise, rng)
        sps = build(Phi, y, "sign", rng); th = sps.theta_ls; W = sps.W; E = ellipsoids(sps)
        if (~E['pd']).sum() >= q: S['skipped'] += 1; continue
        s2 = ((y - Phi @ th) ** 2).sum() / (N - 2); Ri = np.linalg.inv(Phi.T @ Phi); se = np.sqrt(np.diag(Ri) * s2)
        t = time.perf_counter(); Mg, ntg = grid_members(sps, th, se); tg = time.perf_counter() - t
        if len(Mg) == 0: S['skipped'] += 1; continue
        t = time.perf_counter(); pts, nti = radial(sps, th, dirs, W); ti = time.perf_counter() - t
        t = time.perf_counter(); E = ellipsoids(sps); tu = time.perf_counter() - t
        S['t_set']["RACS"] += tg; S['t_set']["IRANY"] += ti; S['t_set']["BUROK"] += tu; S['nt']["RACS"] += ntg; S['nt']["IRANY"] += nti
        thE1, PE1 = rls(Phi, y, R_NOM); thE2, PE2 = rls(Phi, y, s2)
        def hx(Ct):
            Cz = Ct @ W; nrm = np.linalg.norm(Cz, axis=1) + 1e-300
            return support_q(E, Cz / nrm[:, None]) * nrm
        ell = lambda Pm, t_, Sg: np.abs(Pm @ t_ - GDT) + np.sqrt(C95 * np.einsum('ij,jl,il->i', Pm, Sg, Pm))
        fn = {"NONE": None,
              "NOM": lambda Pm: np.abs(Pm @ thn - GDT), "CE": lambda Pm: np.abs(Pm @ th - GDT),
              "ASY": lambda Pm: ell(Pm, th, s2 * Ri), "EKF-nom": lambda Pm: ell(Pm, thE1, PE1), "EKF-becs": lambda Pm: ell(Pm, thE2, PE2),
              "RACS": lambda Pm: np.abs(Pm @ Mg.T - GDT).max(1), "IRANY": lambda Pm: np.abs(Pm @ pts.T - GDT).max(1),
              "BUROK": lambda Pm: np.maximum(Pm @ th + hx(Pm) - GDT, -(Pm @ th) + hx(-Pm) + GDT),
              "ORA": lambda Pm: np.abs(Pm @ th0 - GDT)}
        # lefedes
        d0 = th0 - th; S['cov']["ASY"] += int(d0 @ np.linalg.solve(s2 * Ri, d0) <= C95)
        for nm, (te, Pe) in (("EKF-nom", (thE1, PE1)), ("EKF-becs", (thE2, PE2))):
            dd = th0 - te; S['cov'][nm] += int(dd @ np.linalg.solve(Pe, dd) <= C95)
        S['cov']["RACS"] += int(sps.member(th0[None])[0]); S['cov']["IRANY"] += int(Path(pts).contains_point(th0))
        z0 = C720 @ np.linalg.solve(W, th0 - th); S['cov']["BUROK"] += int((z0 <= support_q(E, C720) + 1e-9).all())
        xs = rng.uniform(-1, 1, K); un = rng.uniform(0.2, 1.8, K)
        for j, nm in enumerate(NAMES):
            if nm == "NONE": u_ = un.copy()
            else:
                t = time.perf_counter(); u_ = robust_u(fn[nm], xs, un); S['t_q'][nm] += time.perf_counter() - t
            nanm = np.isnan(u_); ok = ~nanm
            vv = (np.abs(th0[0] * xs[ok] + th0[1] * u_[ok] - GDT) > V + 1e-12) if ok.any() else np.zeros(0, bool)
            S['viol'][i, j] = vv.mean() if ok.any() else 0.0; S['vany'][i, j] = vv.any() if ok.any() else False
            S['inf'][i, j] = nanm.mean(); S['dev'][j] += list(np.abs(u_[ok] - un[ok]))
        S['nq'] += K; i += 1
    S['n'] = R
    return S


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    conds = sys.argv[2].split(",") if len(sys.argv) > 2 else ["laplace", "szeles", "heavy"]
    Ns = [int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [12, 20, 40]
    print(f"SZIMULACIO (E24, dron fuggoleges tengely). dt={DT} s, V={V} m/s, {R} ablak/cella, {K} dontes/ablak; teher +-30%, legellenallas a=0.95..0.99; EKF prior: nominalis +- (0.02, 0.10); R_nom={R_NOM**0.5:.2f}^2")
    store = {}
    for noise in conds:
        for N in Ns:
            S = run(N, noise, R, np.random.default_rng(30000 + N + 11 * len(noise))); store[(noise, N)] = S
            print(f"\n=== zaj: {noise}, N={N} (T={N*DT:.1f} s) | ablak {R} (kihagyott: vakuum burok {S['skipped']}) | ms/ablak halmaz: racs {S['t_set']['RACS']/R*1e3:.1f}, irany {S['t_set']['IRANY']/R*1e3:.1f}, burok {S['t_set']['BUROK']/R*1e3:.2f} | tesztek/ablak: racs {S['nt']['RACS']/R:.0f}, irany {S['nt']['IRANY']/R:.0f}, burok 0")
            print(f"{'modszer':10s}{'serulo dontes %':>16s}{'van serules ablak % [Wilson]':>32s}{'nincs megold. %':>17s}{'atl |u-unom|':>14s}{'lek. us/dontes':>16s}{'lefedes':>9s}")
            for j, nm in enumerate(NAMES):
                va = S['vany'][:, j].sum(); lo, hi = wilson(int(va), R)
                cv = f"{S['cov'][nm]/R:9.3f}" if nm in S['cov'] else f"{'-':>9s}"
                print(f"{nm:10s}{S['viol'][:, j].mean()*100:16.3f}{va/R*100:14.1f} [{lo*100:.1f},{hi*100:.1f}]{S['inf'][:, j].mean()*100:17.1f}{np.mean(S['dev'][j]) if S['dev'][j] else float('nan'):14.3f}{(S['t_q'][nm]/S['nq']*1e6 if nm != 'NONE' else 0):16.0f}{cv}", flush=True)
    print("\n=== PAROSITOTT TESZT (ablak-szintu 'van serules', McNemar pontos binomialis, ketoldalu p) ===")
    pairs = [("BUROK", "CE"), ("BUROK", "NOM"), ("BUROK", "ASY"), ("BUROK", "EKF-nom"), ("BUROK", "EKF-becs"), ("BUROK", "RACS"), ("IRANY", "RACS"), ("BUROK", "IRANY")]
    for (noise, N), S in store.items():
        print(f"\n{noise}, N={N}:")
        for a, b in pairs:
            ja, jb = NAMES.index(a), NAMES.index(b); x, y_ = S['vany'][:, ja], S['vany'][:, jb]
            b01 = int((x & ~y_).sum()); b10 = int((~x & y_).sum()); p = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 > 0 else 1.0
            print(f"  {a:8s} vs {b:9s}: csak {a} serul {b01:4d}, csak {b} serul {b10:4d} | p = {p:.4f}" + ("  *" if p < 0.05 else ""))
