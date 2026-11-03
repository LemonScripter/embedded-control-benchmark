"""
E28: 500 DRONOS FLOTTA, KOZOS (DRONONKENT NEM HANGOLT) BEALLITAS. SZIMULACIO.
Kerdes: ha a flotta minden dronja UGYANAZT a beallitast kapja (EKF: R_NOM, Q; burok: m=100, q=5, Nw=20, deklarált +-30%),
a zaj pedig dronok kozott/idoben valtozik es nem normalis (vastag farku, zajrohamos), melyik architektura tart tobb dront biztonsagban?
Elore rogzitett forgatokonyvek (a futas elott, nem az eredmeny alapjan valasztva):
  laplace : Laplace zaj, sd 0.05 (a referencia, mint E26)
  vastag  : Student-t(1.5) * 0.03  (vegtelen szoras: rotor-lokes, vibracios tuskek)
  roham   : Laplace sd 0.05, de 3 darab 30 lepeses roham 6x szorassal (nem stacionárius zaj)
  vegyes  : minden dron mas: zajtipus {laplace, t1.5, roham, szeles(0.12)} egyenletesen, lépték exp(N(0,0.4))
Minden dron: teherváltás a 120. lepesben d~U(+-0.3) (ervenyes deklaralt tartomany), uj a. Nw=20. Ugyanaz a szel/dither/referencia minden architekturara (párosítva).
Architekturak: PID, EKF-nom (fix R), EKF-becs (adaptiv R), EKF-rob (adaptiv R + Huber-vagas), BUROK+D, DOBOZ, ORA.
Mertek: sertő lepes % (flotta atlag), 'tartosan serto dron' (>2% sertő lepes) %, a legrosszabb 5% dron atlag serto lepes %, van-serules dron %, beavatkozas, kovetesi hiba.
"""
import sys, time
import numpy as np
from multiprocessing import Pool
from scipy.stats import binomtest
import kolumban_repules as R

T, TCH = R.T, R.TCH
ARCH = ["PID", "EKF-nom", "EKF-becs", "EKF-rob", "BUROK+D", "DOBOZ", "ORA"]
TIPUS = ["laplace", "t15", "roham", "szeles"]


def zaj(tipus, skala, rng):
    if tipus == "laplace": return rng.laplace(size=T) * 0.05 / np.sqrt(2) * skala
    if tipus == "szeles": return rng.laplace(size=T) * 0.12 / np.sqrt(2) * skala
    if tipus == "t15": return rng.standard_t(1.5, size=T) * 0.03 * skala
    w = rng.laplace(size=T) * 0.05 / np.sqrt(2) * skala          # roham
    for _ in range(3):
        s = int(rng.integers(20, T - 30)); w[s:s + 30] *= 6.0
    return w


def epizod(scen, rng):
    ep = R.make_episode("laplace", "ervenyes", rng)               # parameterek, referencia, dither: mint E26
    if scen == "vegyes": tipus = TIPUS[int(rng.integers(4))]; skala = float(np.exp(rng.normal(0, 0.4)))
    else: tipus = {"laplace": "laplace", "vastag": "t15", "roham": "roham"}[scen]; skala = 1.0
    ep["w"] = zaj(tipus, skala, rng)
    return ep


def egy_dron(args):
    scen, seed, Nw = args
    rng = np.random.default_rng(seed); ep = epizod(scen, rng)
    out = {}
    for a in ARCH:
        tm = {k: 0.0 for k in ARCH}
        out[a] = R.run_arch(a, ep, Nw, np.random.default_rng(int(rng.integers(1 << 30))), tm)
    return out


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    Nw = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    scens = sys.argv[3].split(",") if len(sys.argv) > 3 else ["laplace", "vastag", "roham", "vegyes"]
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    print(f"SZIMULACIO (E28, flotta). {N} dron/forgatokonyv, T={T}, teherváltás {TCH}. lepesben, Nw={Nw}, kozos beallitas, ervenyes deklaralt tartomany", flush=True)
    for scen in scens:
        t0 = time.time()
        seeds = [(scen, 810000 + 1000 * len(scen) + i, Nw) for i in range(N)]
        with Pool(nproc) as pool: res = pool.map(egy_dron, seeds, chunksize=4)
        print(f"\n=== forgatokonyv: {scen} | {N} dron | futasido {time.time()-t0:.0f} s")
        print(f"{'architektura':10s}{'serulo lepes %':>15s}{'tartosan serto dron %':>23s}{'legrosszabb 5% dron':>21s}{'van serules dron % [Wilson]':>30s}{'beavatk.':>10s}{'kov. hiba':>10s}{'vesz %':>8s}", flush=True)
        for a in ARCH:
            v = np.array([r[a]['viol'] for r in res]); anyv = sum(r[a]['anyv'] for r in res); lo, hi = wilson(anyv, N)
            tail = np.sort(v)[-max(1, N // 20):].mean()
            print(f"{a:10s}{v.mean()*100:15.3f}{(v > 0.02).mean()*100:23.1f}{tail*100:21.2f}{anyv/N*100:15.1f} [{lo*100:.1f},{hi*100:.1f}]{np.mean([r[a]['dev'] for r in res]):10.3f}{np.mean([r[a]['trk'] for r in res]):10.3f}{np.mean([r[a]['emerg'] for r in res])*100:8.2f}", flush=True)
        print("  PAROSITOTT (dron-szintu 'tartosan serto' >2%, McNemar):")
        for x_, y_ in (("BUROK+D", "EKF-nom"), ("BUROK+D", "EKF-becs"), ("BUROK+D", "EKF-rob"), ("BUROK+D", "DOBOZ"), ("EKF-rob", "EKF-becs"), ("EKF-becs", "PID")):
            a1 = np.array([r[x_]['viol'] > 0.02 for r in res]); a2 = np.array([r[y_]['viol'] > 0.02 for r in res]); b01 = int((a1 & ~a2).sum()); b10 = int((~a1 & a2).sum())
            p = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 else 1.0
            print(f"    {x_:8s} vs {y_:9s}: csak {x_} serto {b01:4d}, csak {y_} serto {b10:4d} | p = {p:.4f}" + ("  *" if p < 0.05 else ""), flush=True)
