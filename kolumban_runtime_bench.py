"""
Kolumban-tetel (SPS, regresszios eset) futasideju koltsege - E3 meres.
Kerdes: hany ms a) naiv, b) jol megirt SPS, c) racsos regio-feltaras,
        d) iranyitott (radialis bisection) feltaras, e) futasideju ellenorzes.
Modell: y_t = phi_t' theta0 + n_t  (meresi allapotu regresszio, x+ = Theta*phi + w)
Szint: 1 - q/m.  Merje a PC-n a ms-ot; ESP32-re csak MAC-szambol becsulunk.
"""
import time
import numpy as np

rng = np.random.default_rng(1)


def make_problem(n_theta=3, N=50, m=100, q=5, noise=0.3):
    theta0 = rng.normal(size=n_theta)
    Phi = rng.normal(size=(N, n_theta))               # exogen regresszor
    y = Phi @ theta0 + noise * rng.laplace(size=N)    # szimmetrikus, nem gauss zaj
    alpha = rng.choice([-1.0, 1.0], size=(m - 1, N))  # m-1 elojelsorozat
    return theta0, Phi, y, alpha, m, q


class SPS:
    """Hatekony SPS: elofeldolgozas ablakonkent, theta-nkent O(m n^2)."""

    def __init__(self, Phi, y, alpha, m, q):
        self.m, self.q = m, q
        N, n = Phi.shape
        R = Phi.T @ Phi / N
        w, V = np.linalg.eigh(R)
        self.W = V @ np.diag(w ** -0.5) @ V.T        # R^{-1/2}, egyszer
        # g0(theta) = b0 - R theta ; gi(theta) = bi - Ai theta
        self.b0 = Phi.T @ y / N
        self.A0 = R
        self.B = (alpha * y) @ Phi / N               # (m-1, n)
        self.A = np.einsum('it,tj,tk->ijk', alpha, Phi, Phi) / N  # (m-1,n,n)
        self.theta_ls = np.linalg.solve(R, self.b0)

    def member(self, thetas):
        """thetas: (P, n) -> bool (P,). Batch, vektorizalt."""
        th = np.atleast_2d(thetas)
        g0 = (self.b0 - th @ self.A0.T) @ self.W                   # (P,n)
        gi = self.B[None] - np.einsum('ijk,pk->pij', self.A, th)   # (P,m-1,n)
        gi = gi @ self.W
        S0 = (g0 ** 2).sum(-1)
        Si = (gi ** 2).sum(-1)
        n_bigger = (Si > S0[:, None]).sum(-1)       # hany perturbalt nagyobb
        return n_bigger >= self.q                   # S0 nincs a q legnagyobb kozott


def naive_member(Phi, y, alpha, q, theta):
    """Demo-stilus: perturbaciononkent ujra inv(), ciklusban."""
    N, n = Phi.shape
    def S(a):
        eps = y - Phi @ theta
        g = Phi.T @ (a * eps)
        Rinv = np.linalg.inv(Phi.T @ Phi)           # minden korben ujra (strawman)
        return g @ Rinv @ g
    S0 = S(np.ones(N))
    return sum(S(a) > S0 for a in alpha) >= q


def timeit(f, reps=20):
    f()
    t = time.perf_counter()
    for _ in range(reps):
        f()
    return (time.perf_counter() - t) / reps * 1e3


def radial_extent(sps, center, u, r0=0.05, iters=22):
    """Csillag-konvex regio feltetelezesevel: legnagyobb r, hogy center+r*u tagja."""
    r = r0
    tests = 0
    while sps.member(center + r * u)[0] and r < 1e3:
        r *= 2
        tests += 1
    lo, hi = r / 2, r
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        tests += 1
        if sps.member(center + mid * u)[0]:
            lo = mid
        else:
            hi = mid
    return lo, tests + 1


def run(n_theta, m=100, q=5, N=50, grid_pts=41):
    theta0, Phi, y, alpha, m, q = make_problem(n_theta, N, m, q)
    t_pre = timeit(lambda: SPS(Phi, y, alpha, m, q), reps=5)
    sps = SPS(Phi, y, alpha, m, q)
    c = sps.theta_ls

    print(f"\n=== n_theta={n_theta}, N={N}, m={m}, q={q}  (szint {1-q/m:.2f}) ===")
    t_naive = timeit(lambda: naive_member(Phi, y, alpha, q, theta0), reps=3)
    t_one = timeit(lambda: sps.member(theta0), reps=200)
    macs_one = (m - 1) * n_theta * n_theta + (m - 1) * n_theta * n_theta  # A.th + .W
    print(f"naiv SPS, 1 theta-teszt (inv ciklusban) : {t_naive:9.3f} ms")
    print(f"hatekony SPS, 1 theta-teszt (vektoriz.) : {t_one:9.3f} ms   (~{macs_one} MAC)")
    print(f"ablak-elofeldolgozas (A_i, b_i, R^-1/2) : {t_pre:9.3f} ms   (ablakonkent, nem lepesenkent)")

    # racsos feltaras: a teszt-szamot szamoljuk; ha tul sok, extrapolalunk
    P = grid_pts ** n_theta
    spread = 4 * np.sqrt(np.diag(np.linalg.inv(Phi.T @ Phi / N)) / N) * 1.0 + 0.5
    if P <= 300_000:
        axes = [np.linspace(c[j] - spread[j] * 3, c[j] + spread[j] * 3, grid_pts) for j in range(n_theta)]
        G = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, n_theta)
        t0 = time.perf_counter()
        ok = np.concatenate([sps.member(G[i:i + 20000]) for i in range(0, len(G), 20000)])
        t_grid = (time.perf_counter() - t0) * 1e3
        print(f"racs ({grid_pts}^{n_theta}={P} teszt)            : {t_grid:9.1f} ms   (elfogadott: {ok.sum()})")
        true_max = np.array([G[ok][:, j].max() for j in range(n_theta)])
        true_min = np.array([G[ok][:, j].min() for j in range(n_theta)])
    else:
        t_grid = t_one * P
        print(f"racs ({grid_pts}^{n_theta}={P:.2e} teszt)       : {t_grid:9.1f} ms   (EXTRAPOLALT, nem futtatva)")
        true_max = true_min = None

    # iranyitott: +-tengely iranyok (K=2n)
    t0 = time.perf_counter()
    ext, tests = [], 0
    for j in range(n_theta):
        for s in (+1, -1):
            u = np.zeros(n_theta); u[j] = s
            r, t = radial_extent(sps, c, u)
            ext.append((j, s, r)); tests += t
    t_dir = (time.perf_counter() - t0) * 1e3
    print(f"iranyitott (2n={2*n_theta} irany, {tests} teszt)    : {t_dir:9.3f} ms")
    if true_max is not None:
        gap = []
        for j, s, r in ext:
            tm = (true_max[j] - c[j]) if s > 0 else (c[j] - true_min[j])
            gap.append(r - tm)
        print(f"  radialis - racs-kiterjedes (negativ = a bisection ALULBECSUL): min={min(gap):+.3f} max={max(gap):+.3f}")
        print(f"  (racs-felbontas: {2*3*spread.max()/(grid_pts-1):.3f})")
    box_lo = np.array([c[j] - [r for jj, s, r in ext if jj == j and s < 0][0] for j in range(n_theta)])
    box_hi = np.array([c[j] + [r for jj, s, r in ext if jj == j and s > 0][0] for j in range(n_theta)])
    print(f"  theta0 a radialis dobozban: {bool(np.all(theta0 >= box_lo) and np.all(theta0 <= box_hi))}")

    # futasideju ellenorzes: legrosszabb c'theta az ellipszoid-burkon (n^2 MAC)
    Sig = np.linalg.inv(Phi.T @ Phi) * (np.mean((y - Phi @ c) ** 2)) * 6.0
    L = np.linalg.cholesky(Sig)
    cvec = rng.normal(size=n_theta)
    def runtime_check():
        return cvec @ c - np.linalg.norm(L.T @ cvec)
    t_rt = timeit(runtime_check, reps=2000) * 1e3
    print(f"futasideju ellenorzes (ellipszoid, kenyszerenkent): {t_rt:9.2f} us   (~{n_theta*n_theta+n_theta} MAC)")

    # ESP32 becsles (FELTEVES: ~10 MFLOP/s tartosan, f32 FPU, 2 FLOP/MAC)
    mfl = 10e6
    def esp(macs): return macs * 2 / mfl * 1e3
    print(f"ESP32 BECSLES (feltevesen 10 MFLOP/s): 1 theta-teszt ~{esp(macs_one):.2f} ms, "
          f"iranyitott feltaras ~{esp(macs_one*tests):.1f} ms, futasideju ellenorzes ~{esp(n_theta*n_theta+n_theta)*1e3:.1f} us")
    return dict(naive=t_naive, one=t_one, grid=t_grid, direct=t_dir, rt_us=t_rt, tests=tests)


if __name__ == '__main__':
    for n in (3, 6, 12):
        run(n)
