"""
ARX-eset: y_t = a1 y_{t-1} + a2 y_{t-2} + b u_{t-1} + w_t.  A regresszor a MULTBELI KIMENETEKET tartalmazza
-> a 2.17. tetel 'exogen regresszor' feltetele serul. Merjuk:
 (A) a pontos SPS/DP-teszt lefedettsege (theta0 a regioban?)   [elvart 0.95]
 (B) az L3 tartalmazas (algebrai, nem fugg a zajtol): theta0 a burokban?  (szukseges feltetel-ellenorzes)
 (C) a burok vakuum-aranya, szorossag
"""
import sys
import numpy as np
from kolumban_runtime_bench import SPS

q, m = 5, 100


def gen(rng, a, b, N, noise, burn=200):
    T = N + burn + 2
    u = rng.normal(size=T)
    if noise == "gauss": w = rng.normal(size=T) * 0.5
    elif noise == "laplace": w = rng.laplace(size=T) * 0.35
    else: w = (rng.exponential(size=T) - 1.0) * 0.5          # ferde, nem-zerus-szimmetrikus (de zerus atlag)
    y = np.zeros(T)
    for t in range(2, T):
        y[t] = a[0] * y[t-1] + a[1] * y[t-2] + b * u[t-1] + w[t]
    Phi = np.stack([y[1:-1], y[:-2], u[1:-1]], 1)[burn:burn+N]
    yy = y[2:][burn:burn+N]
    return Phi, yy


def build(Phi, y, mode, rng):
    N = len(y)
    al = rng.choice([-1.0, 1.0], size=(m - 1, N))
    sps = SPS(Phi, y, al, m, q)
    if mode == "perm":
        P = np.array([rng.permutation(N) for _ in range(m - 1)])
        sps.B = np.einsum('tj,it->ij', Phi, y[P]) / N
        sps.A = np.einsum('tj,itk->ijk', Phi, Phi[P]) / N
    return sps


def ellipsoids(sps):
    W = sps.W; th = sps.theta_ls; n = len(th)
    e = sps.B - np.einsum('ijk,k->ij', sps.A, th)
    f = e @ W.T
    M = np.einsum('ab,ibc,cd->iad', W, sps.A, W)
    Pm = np.eye(n)[None] - np.einsum('iba,ibc->iac', M, M)
    lam = np.linalg.eigvalsh(Pm)[:, 0]
    pd = lam > 1e-9
    Pinv = np.zeros_like(Pm); Pinv[pd] = np.linalg.inv(Pm[pd])
    u = np.einsum('iba,ib->ia', M, f)
    zeta = -np.einsum('iab,ib->ia', Pinv, u)
    kappa = (f ** 2).sum(1) + np.einsum('ia,iab,ib->i', u, Pinv, u)
    return dict(zeta=zeta, Pinv=Pinv, kappa=kappa, pd=pd)


def support_q(E, C):
    t1 = C @ np.where(E['pd'][:, None], E['zeta'], 0).T
    quad = np.einsum('kj,ijl,kl->ki', C, E['Pinv'], C)
    h = t1 + np.sqrt(np.maximum(np.where(E['pd'], E['kappa'], 0), 0)[None] * np.maximum(quad, 0))
    h = np.where(E['pd'][None], h, np.inf)
    return np.sort(h, axis=1)[:, ::-1][:, q - 1]


def run(a, b, N, noise, mode, R, rng):
    th0 = np.array([a[0], a[1], b])
    ex = hull = vac = 0; n_ok = 0
    ratios = []
    for _ in range(R):
        Phi, y = gen(rng, a, b, N, noise)
        sps = build(Phi, y, mode, rng)
        ex += int(sps.member(th0[None])[0])
        E = ellipsoids(sps)
        if (~E['pd']).sum() >= q:
            vac += 1; hull += 1; continue
        z0 = np.linalg.solve(sps.W, th0 - sps.theta_ls)
        C = rng.normal(size=(60, 3)); C = np.vstack([C, z0[None]])
        C /= np.linalg.norm(C, axis=1)[:, None]
        h = support_q(E, C)
        hull += int(np.all(C @ z0 <= h + 1e-9))
        ratios.append(np.median(h / np.maximum(C @ z0, 1e-9)))
    return ex / R, hull / R, vac / R


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    cases = {"a=(.5,-.2) mild": (0.5, -0.2), "a=(.9,-.1) near-unit": (0.9, -0.1), "a=(1.2,-.5) resonant": (1.2, -0.5)}
    print(f"R={R}, m={m}, q={q}  cel: pontos teszt ~0.95 (0.95 = 1-q/m)")
    print(f"{'rendszer':24s}{'zaj':8s}{'N':>5s}{'mod':>6s} | {'pontos teszt':>12s} {'burok':>7s} {'vakuum':>7s}")
    for name, a in cases.items():
        for noise in ("gauss", "laplace", "skew"):
            for N in (50, 200):
                for mode in ("sign", "perm"):
                    c, h, v = run(a, 1.0, N, noise, mode, R, rng)
                    print(f"{name:24s}{noise:8s}{N:5d}{mode:>6s} | {c:12.3f} {h:7.3f} {v:7.2f}", flush=True)
