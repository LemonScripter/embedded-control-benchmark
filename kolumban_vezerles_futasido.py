"""
E23: A KOLUMBAN-TETEL A VEZERLESBEN, INFORMACIOSZEGENY KORNYEZETBEN (N=12..40): mibe kerul a garancia, es az irany + lapitas futasidore csokkenti-e? SZIMULACIO.
Ugyanaz a feladat, mint az E10/E12/E13: instabil x+ = a x + b u + w (a0=1.1, b0=1.0), theta=(a,b) ismeretlen, biztonsagi szuro |x+| <= 1.5 MINDEN hihetonek tartott theta-ra.
Ugyanabbol az ablakbol, parositva, harom modja a "hihetoek halmazanak" eloallitasara + a szuro lekerdezese:
  RACS   : a pontos SPS-halmaz racson (121 x 121 = 14 641 teszt, ha a halmaz kilog a racsbol: a racs duplazodik) -> szuro: a halmaz pontjain a legrosszabb eset  (a tetel "drag" utja)
  IRANY  : 16 irany, iranyonkenti duplazas + 14 felezes -> a hatarpontok sokszoge (belso kozelites, nincs garancia)
  BUROK  : lapitas, zart alak, 0 teszt (kulso kozelites, tartalmazza a halmazt)
Mertek: (1) ido/ablak: az ablak SPS-elofeldolgozasa kozos (m-1 elojelsorozat); a halmaz-eloallitas ideje; (2) szuro-lekerdezes ideje (us/dontes: 161 jeloltre a legrosszabb eset);
        (3) a VALODI theta0 lefedese a halmazzal (cel 0.95; zart hurkon nem bizonyitott), (4) serulo dontes %, nincs megoldas %, atl |u-unom| (E13 mertekei).
"""
import sys, time
import numpy as np
from matplotlib.path import Path
from kolumban_arx import build, ellipsoids, support_q, q, m
import kolumban_decision as kd

a0, b0, X = kd.a0, kd.b0, kd.X
ND = 16; NG = 121; K = 40


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def grid_members(sps, th, se):
    """pontos SPS-halmaz racson; a racs duplazodik, amig a halmaz nem er a racs szelehez. -> (tagok, tesztek szama)"""
    H = 8 * se; nt = 0
    for _ in range(6):
        gx = np.linspace(-1, 1, NG)
        GA, GB = np.meshgrid(gx * H[0], gx * H[1]); P = th + np.stack([GA.ravel(), GB.ravel()], 1)
        mem = np.concatenate([sps.member(P[s:s + 4000]) for s in range(0, len(P), 4000)]); nt += len(P)
        M2 = mem.reshape(NG, NG)
        if not (M2[[0, -1], :].any() or M2[:, [0, -1]].any()): return P[mem], nt
        H = H * 2
    return P[mem], nt


def radial(sps, th, dirs, W):
    nt = 0; pts = []
    for d in dirs:
        dd = W @ d; r = 0.02
        while sps.member((th + r * dd)[None])[0] and r < 1e3: r *= 2; nt += 1
        lo, hi = r / 2, r
        for _ in range(14):
            mid = 0.5 * (lo + hi); nt += 1
            if sps.member((th + mid * dd)[None])[0]: lo = mid
            else: hi = mid
        pts.append(th + lo * dd)
    return np.array(pts), nt


def run(N, noise, R, rng):
    kd.rng = rng
    ph = np.linspace(0, 2 * np.pi, ND, endpoint=False); dirs = np.stack([np.cos(ph), np.sin(ph)], 1)
    ang = np.linspace(0, 2 * np.pi, 720, endpoint=False); C720 = np.stack([np.cos(ang), np.sin(ang)], 1)
    names = ["RACS", "IRANY", "BUROK"]
    S = dict(n=0, t_build=0.0, t_set={k: 0.0 for k in names}, t_q={k: 0.0 for k in names}, nt={"RACS": [], "IRANY": []},
             cov={k: 0 for k in names}, vf={k: 0.0 for k in names}, inf={k: 0.0 for k in names}, dev={k: [] for k in names}, nq=0)
    th0 = np.array([a0, b0])
    for _ in range(R):
        Phi, y = kd.window(N, noise)
        t = time.perf_counter(); sps = build(Phi, y, "sign", rng); tb = time.perf_counter() - t
        th = sps.theta_ls; W = sps.W
        E = ellipsoids(sps)
        if (~E['pd']).sum() >= q: continue
        s2 = ((y - Phi @ th) ** 2).sum() / (N - 2); se = np.sqrt(np.diag(np.linalg.inv(Phi.T @ Phi)) * s2)
        t = time.perf_counter(); Mg, ntg = grid_members(sps, th, se); tg = time.perf_counter() - t
        if len(Mg) == 0: continue
        t = time.perf_counter(); pts, nti = radial(sps, th, dirs, W); ti = time.perf_counter() - t
        t = time.perf_counter(); E = ellipsoids(sps); tu = time.perf_counter() - t
        def hull_extent(Ct):
            Cz = Ct @ W; nrm = np.linalg.norm(Cz, axis=1) + 1e-300
            return support_q(E, Cz / nrm[:, None]) * nrm
        fn = {"RACS": lambda Pm: np.abs(Pm @ Mg.T).max(1),
              "IRANY": lambda Pm: np.abs(Pm @ pts.T).max(1),
              "BUROK": lambda Pm: np.maximum(Pm @ th + hull_extent(Pm), -(Pm @ th) + hull_extent(-Pm))}
        S['n'] += 1; S['t_build'] += tb; S['t_set']["RACS"] += tg; S['t_set']["IRANY"] += ti; S['t_set']["BUROK"] += tu
        S['nt']["RACS"].append(ntg); S['nt']["IRANY"].append(nti)
        # lefedes: a valodi theta0 a halmazban van-e
        S['cov']["RACS"] += int(sps.member(th0[None])[0])
        S['cov']["IRANY"] += int(Path(pts).contains_point(th0))
        z0 = C720 @ np.linalg.solve(W, th0 - th); S['cov']["BUROK"] += int((z0 <= support_q(E, C720) + 1e-9).all())
        xs = rng.uniform(-1, 1, K); un = rng.uniform(-2.5, 2.5, K)
        for nm in names:
            t = time.perf_counter(); u_ = kd.robust_u(None, fn[nm], xs, un); S['t_q'][nm] += time.perf_counter() - t
            nanm = np.isnan(u_); ok = ~nanm
            S['inf'][nm] += nanm.mean()
            S['vf'][nm] += (np.abs(a0 * xs[ok] + b0 * u_[ok]) > X + 1e-12).mean() if ok.any() else 0.0
            S['dev'][nm] += list(np.abs(u_[ok] - un[ok]))
        S['nq'] += K
    return S


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    conds = sys.argv[2].split(",") if len(sys.argv) > 2 else ["laplace", "hetero", "heavy"]
    Ns = [int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [12, 20, 40]
    print(f"SZIMULACIO (E23). a0={a0} b0={b0} X={X}; {R} ablak/cella; racs {NG}x{NG}; irany {ND}. Cel lefedes 0.95 (zart hurkon nem bizonyitott).")
    store = {}
    for noise in conds:
        for N in Ns:
            S = run(N, noise, R, np.random.default_rng(9000 + N + 11 * len(noise))); store[(noise, N)] = S; n = max(S['n'], 1)
            print(f"\n=== zaj: {noise}, N={N} | ervenyes ablak {S['n']}/{R} | kozos SPS-elofeldolgozas {S['t_build']/n*1e3:.2f} ms | tesztek/ablak: racs {np.mean(S['nt']['RACS']):.0f}, irany {np.mean(S['nt']['IRANY']):.0f}, burok 0")
            print(f"{'halmaz':8s}{'ido/ablak ms':>14s}{'lekerdezes us/dontes':>22s}{'lefedes':>18s}{'serulo %':>10s}{'nincs megold. %':>17s}{'atl |u-unom|':>14s}")
            for nm in ("RACS", "IRANY", "BUROK"):
                c = S['cov'][nm]; lo, hi = wilson(c, n)
                print(f"{nm:8s}{S['t_set'][nm]/n*1e3:14.2f}{S['t_q'][nm]/max(S['nq'],1)*1e6:22.0f}{c/n:9.3f} [{lo:.3f},{hi:.3f}]{S['vf'][nm]/n*100:10.2f}{S['inf'][nm]/n*100:17.1f}{np.mean(S['dev'][nm]) if S['dev'][nm] else float('nan'):14.3f}", flush=True)
    print("\n=== OSSZEVONT zajtipusokra, N szerint ===")
    for N in Ns:
        ks = [(nz, N) for nz in conds]; n = sum(store[k]['n'] for k in ks)
        print(f"\nN={N} (n={n}): ido/ablak (ms): " + ", ".join(f"{nm} {sum(store[k]['t_set'][nm] for k in ks)/n*1e3:.2f}" for nm in ("RACS", "IRANY", "BUROK"))
              + f" | kozos elofeldolgozas {sum(store[k]['t_build'] for k in ks)/n*1e3:.2f}"
              + " | lekerdezes us/dontes: " + ", ".join(f"{nm} {sum(store[k]['t_q'][nm] for k in ks)/sum(store[k]['nq'] for k in ks)*1e6:.0f}" for nm in ("RACS", "IRANY", "BUROK"))
              + " | lefedes: " + ", ".join(f"{nm} {sum(store[k]['cov'][nm] for k in ks)/n:.3f}" for nm in ("RACS", "IRANY", "BUROK"))
              + " | serulo %: " + ", ".join(f"{nm} {sum(store[k]['vf'][nm] for k in ks)/n*100:.2f}" for nm in ("RACS", "IRANY", "BUROK"))
              + " | nincs megoldas %: " + ", ".join(f"{nm} {sum(store[k]['inf'][nm] for k in ks)/n*100:.1f}" for nm in ("RACS", "IRANY", "BUROK")))
