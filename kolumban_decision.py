"""
DONTES-teszt: javitja-e a Kolumban-burok a dontest informaciohianyos kornyezetben?
Plant: x+ = a0 x + b0 u + w  (a0=1.1 instabil), zart hurok az adatgyujtesre (u=-k x+d), N mintas ablak, theta=(a,b) ismeretlen.
Biztonsagi szuro: a legkozelebbi u az u_nom-hoz, amire |a x + b u| <= X MINDEN theta-ra a bizonytalansagi halmazban.
Dontesi hiba = a VALODI theta0-ra serul a korlat. Osszehasonlitas:
  CE   : csak theta^ (certainty equivalence, bizonytalansag nelkul)
  ASY  : aszimptotikus 95%-os ellipszoid (klasszikus, normalitasra epul)
  KOL  : Kolumban-L3 burok (elojeles SPS, q-adik legnagyobb ellipszoid-tartoszam), 1-q/m=0.95
  ORA  : a valodi theta0 (nem megvalosithato; a konzervativizmus viszonyitasa)
"""
import sys
import numpy as np
from scipy.stats import chi2
from kolumban_arx import build, ellipsoids, support_q, q, m

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
a0, b0, k_, sd_d, X = 1.1, 1.0, 1.0, 0.5, 1.5
UG = np.linspace(-4, 4, 161)


def window(N, noise):
    sd = np.ones(N)*0.5
    if noise == "hetero": sd = 0.1*np.exp(np.linspace(0, np.log(13), N))      # 13:1 szorasprofil, szimmetrikus
    w = (rng.standard_t(1.5, size=N)*0.3) if noise == "heavy" else rng.laplace(size=N)*sd/np.sqrt(2)
    d = rng.normal(size=N)*sd_d
    x = np.zeros(N+1); u = np.zeros(N)
    for t in range(N):
        u[t] = -k_*x[t] + d[t]; x[t+1] = a0*x[t] + b0*u[t] + w[t]
    return np.stack([x[:-1], u], 1), x[1:]


def robust_u(phi_fn, worst, xs, unom):
    """minden (x,unom) parra: a legkozelebbi UG-beli u, amire worst(x,u) <= X; nan ha nincs"""
    out = np.full(len(xs), np.nan)
    for i, (x, un) in enumerate(zip(xs, unom)):
        Phi = np.stack([np.full_like(UG, x), UG], 1)
        ok = worst(Phi) <= X
        if ok.any():
            j = np.argmin(np.where(ok, np.abs(UG-un), np.inf)); out[i] = UG[j]
    return out


def run(N, noise, R, K=40):
    keys = ("CE", "ASY", "KOL", "ORA", "NONE")
    viol_any = {k: 0 for k in keys}; viol_frac = {k: 0.0 for k in keys}
    infeas = {k: 0.0 for k in keys}; dev = {k: [] for k in keys}; vac = 0
    for _ in range(R):
        Phi, y = window(N, noise)
        sps = build(Phi, y, "sign", rng); th = sps.theta_ls; W = sps.W
        E = ellipsoids(sps)
        vacuous = (~E['pd']).sum() >= q; vac += vacuous
        rss = ((y-Phi@th)**2).sum(); s2 = rss/(N-2); Ri = np.linalg.inv(Phi.T@Phi)
        c95 = chi2.ppf(0.95, 2)
        xs = rng.uniform(-1, 1, K); un = rng.uniform(-2.5, 2.5, K)

        def w_ce(P): return np.abs(P@th)
        def w_asy(P): return np.abs(P@th) + np.sqrt(c95*s2*np.einsum('ij,jl,il->i', P, Ri, P))
        def w_kol(P):
            # max_{theta in hull} +-theta'phi = +-th'phi + h_q(W' (+-phi))   (z-ter)
            out = np.empty(len(P))
            C = P@W                                # (n_u, 2): c = W' phi (W szimmetrikus)
            nrm = np.linalg.norm(C, axis=1)+1e-300
            Cu = C/nrm[:, None]
            hp = support_q(E, Cu)*nrm; hm = support_q(E, -Cu)*nrm
            return np.maximum(P@th + hp, -(P@th) + hm)
        def w_ora(P): return np.abs(P@np.array([a0, b0]))
        cases = {"CE": w_ce, "ASY": w_asy, "ORA": w_ora}
        if not vacuous: cases["KOL"] = w_kol
        sel = {}
        for name, fn in cases.items():
            sel[name] = robust_u(None, fn, xs, un)
        sel["NONE"] = un.copy()
        for name in keys:
            if name not in sel:                     # KOL vakuum: nincs garancia -> infeasible mindenre
                infeas[name] += 1.0; continue
            u_ = sel[name]; nanm = np.isnan(u_)
            infeas[name] += nanm.mean()
            ok = ~nanm
            v = np.abs(a0*xs[ok] + b0*u_[ok]) > X + 1e-12
            viol_frac[name] += v.mean() if ok.any() else 0.0
            viol_any[name] += int(v.any())
            dev[name] += list(np.abs(u_[ok]-un[ok]))
    return {k: (viol_any[k]/R, viol_frac[k]/R, infeas[k]/R, np.mean(dev[k]) if dev[k] else np.nan) for k in keys}, vac/R


if __name__ == "__main__":
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    print(f"plant a0={a0} b0={b0} (nyiltan instabil), X={X}; R={R}; cel: a szuro utan a VALODI plantra serul-e |a x+b u|<=X")
    print("szuro: CE=csak theta^, ASY=aszimptotikus 95%, KOL=Kolumban-burok 95%, ORA=igazi theta0, NONE=nincs szuro")
    print(f"{'N':>4s} {'zaj':>7s} | {'mod':>4s} {'trial-ben van serules':>22s} {'serulo dontesek %':>18s} {'nincs megoldas %':>17s} {'atl |u-unom|':>13s}")
    for noise in (sys.argv[3].split(",") if len(sys.argv) > 3 else ("laplace", "hetero")):
        for N in (12, 20, 40):
            res, vac = run(N, noise, R)
            for name in ("NONE", "CE", "ASY", "KOL", "ORA"):
                va, vf, inf, d = res[name]
                print(f"{N:4d} {noise:>7s} | {name:>4s} {va:22.3f} {vf*100:18.2f} {inf*100:17.1f} {d:13.3f}" + (f"   (KOL vakuum: {vac*100:.0f}%)" if name == "KOL" else ""), flush=True)
