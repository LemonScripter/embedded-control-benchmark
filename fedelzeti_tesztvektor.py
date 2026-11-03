"""Tesztvektor a fedelzeti (Rust, f32) megvalositas ellenorzesehez: egy valodi szimulalt dron-ablak + elojelsorozatok + a Python (f64) referencia-eredmenyek.
Kimenet: fedelzeti/tv.txt"""
import os
import numpy as np
from kolumban_arx import ellipsoids, support_q, q, m
from kolumban_runtime_bench import SPS
import kolumban_dron as D
import kolumban_dron_elozetes as E

N = int(os.environ.get("NW", "20"))
rng = np.random.default_rng(2026)
while True:
    Phi, y, th0 = E.window(N, "laplace", "ervenyes", rng)
    al = rng.choice([-1.0, 1.0], size=(m - 1, N)); sps = SPS(Phi, y, al, m, q); Eel = ellipsoids(sps)
    if (~Eel['pd']).sum() == 0: break
th, W = sps.theta_ls, sps.W


def hx(Ct):
    Cz = Ct @ W; nrm = np.linalg.norm(Cz, axis=1) + 1e-300
    return support_q(Eel, Cz / nrm[:, None]) * nrm


GDT, V, UG = D.GDT, D.V, np.linspace(0, 3, 121)
lines = [f"N {N}", "S"] + [f"{Phi[t,0]:.9g} {Phi[t,1]:.9g} {y[t]:.9g}" for t in range(N)] + ["G"] + [" ".join(str(int(a)) for a in al[i]) for i in range(m - 1)]
lines += [f"T {th[0]:.9g} {th[1]:.9g}", f"W {W[0,0]:.9g} {W[0,1]:.9g} {W[1,1]:.9g}", f"PD {int((~Eel['pd']).sum())}"]
ang = np.linspace(0, 2 * np.pi, 16, endpoint=False); C = np.stack([np.cos(ang), np.sin(ang)], 1); h = hx(C)
lines.append(f"H {len(C)}"); lines += [f"{C[k,0]:.9g} {C[k,1]:.9g} {h[k]:.9g}" for k in range(len(C))]
gate = []
for v in (-0.9, -0.3, 0.3, 0.8, 0.95):
    for un in (0.2, 0.9, 1.6, 2.4):
        Pm = np.stack([np.full_like(UG, v), UG], 1)
        worst = np.maximum(Pm @ th + hx(Pm) - GDT, -(Pm @ th) + hx(-Pm) + GDT); ok = worst <= V
        ue = float(UG[np.argmin(np.where(ok, np.abs(UG - un), np.inf))]) if ok.any() else float('nan')
        gate.append((v, un, ue))
lines.append(f"Q {len(gate)}"); lines += [f"{v:.9g} {un:.9g} {ue:.9g}" for v, un, ue in gate]
os.makedirs("fedelzeti", exist_ok=True)
open(os.path.join("fedelzeti", f"tv{N}.txt"), "w").write("\n".join(lines) + "\n")
print("tv kesz: N =", N, "theta", th, "| a valodi theta", th0, "| gate esetek:", len(gate), "ebbol vegtelen:", sum(1 for g in gate if g[2] != g[2]))
