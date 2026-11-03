"""
E27b: ESP32-C3 (RISC-V RV32IMC, 160 MHz, FPU NELKUL) cel: a Kolumban-burok kapu EMULATORBAN futtatva (Unicorn), UTASITASSZAMLALASSAL.
MERT: az utasitasok szama (tipusonkent), a kod merete, a statikus RAM. A valodi tesztvektoron (ugyanaz, mint a PC-s ellenorzes) a kimenet egyezik a Python (f64) referenciaval.
BECSULT: az ido (ciklus/utasitas feltevesekkel; lent). NEM hardveres meres. A (nem emulalhato) Xtensa ESP32 (LX6, van FPU-ja) ettol kulon, becsult.
Hasznalat: python fedelzeti_c3_meres.py N   (N = 12 | 20 | 40, a fedelzeti/tvN.txt alapjan)
"""
import struct, sys, math
from elftools.elf.elffile import ELFFile
from unicorn import Uc, UC_ARCH_RISCV, UC_MODE_RISCV32, UC_HOOK_CODE
from unicorn.riscv_const import UC_RISCV_REG_A0, UC_RISCV_REG_A1, UC_RISCV_REG_A2, UC_RISCV_REG_SP, UC_RISCV_REG_RA

ELF = "fedelzeti_c3/target/riscv32imc-unknown-none-elf/release/kolumban_c3"
FLASH, FLASH_SZ = 0x42000000, 4 * 1024 * 1024
RAM, RAM_SZ = 0x3FC80000, 384 * 1024
HALT = FLASH + FLASH_SZ - 0x1000

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
NMAX, M1 = 40, 99


def read_tv(path):
    tok = open(path, encoding="utf-8-sig").read().split(); i = 0
    def nx():
        nonlocal i; i += 1; return tok[i - 1]
    assert nx() == "N"; n = int(nx()); assert nx() == "S"
    samples = [[float(nx()), float(nx()), float(nx())] for _ in range(n)]
    assert nx() == "G"; sign = [[int(nx()) for _ in range(n)] for _ in range(M1)]
    assert nx() == "T"; th = [float(nx()), float(nx())]
    assert nx() == "W"; w = [float(nx()), float(nx()), float(nx())]
    assert nx() == "PD"; pd = int(nx())
    assert nx() == "H"; k = int(nx()); hx = [[float(nx()), float(nx()), float(nx())] for _ in range(k)]
    assert nx() == "Q"; k = int(nx()); gate = [[float(nx()), float(nx()), float(nx())] for _ in range(k)]
    return dict(n=n, samples=samples, sign=sign, th=th, w=w, pd=pd, hx=hx, gate=gate)


f2i = lambda x: struct.unpack("<I", struct.pack("<f", x))[0]
i2f = lambda b: struct.unpack("<f", struct.pack("<I", b & 0xFFFFFFFF))[0]

elf = ELFFile(open(ELF, "rb"))
sym = {s.name: s["st_value"] for s in elf.get_section_by_name(".symtab").iter_symbols() if s.name}
img = bytearray(FLASH_SZ)
uc = Uc(UC_ARCH_RISCV, UC_MODE_RISCV32)
uc.mem_map(FLASH, FLASH_SZ); uc.mem_map(RAM, RAM_SZ)
for seg in elf.iter_segments():
    if seg["p_type"] != "PT_LOAD" or seg["p_filesz"] == 0: continue
    data = seg.data(); va = seg["p_vaddr"]
    uc.mem_write(va, data)
    if FLASH <= va < FLASH + FLASH_SZ: img[va - FLASH:va - FLASH + len(data)] = data
sections = {s.name: s for s in elf.iter_sections()}
code_size = sections[".text"]["sh_size"]; ro_size = sections[".rodata"]["sh_size"] if ".rodata" in sections else 0
bss_size = sections[".bss"]["sh_size"] if ".bss" in sections else 0; data_size = sections[".data"]["sh_size"] if ".data" in sections else 0
uc.mem_write(HALT, b"\x01\x00" * 4)   # nop-ok; a HALT cimre erkezve az emulacio leall


def classify(addr):
    off = addr - FLASH; h = img[off] | (img[off + 1] << 8)
    if (h & 3) != 3:        # tomoritett (16 bites)
        quad, f3 = h & 3, h >> 13
        if quad == 0: return ("load" if f3 == 2 else "store" if f3 == 6 else "alu"), 2
        if quad == 1:
            if f3 in (1, 5): return "jump", 2
            if f3 in (6, 7): return "branch", 2
            return "alu", 2
        if f3 == 2: return "load", 2
        if f3 == 6: return "store", 2
        if f3 == 4:
            rs2 = (h >> 2) & 31; rs1 = (h >> 7) & 31; b12 = (h >> 12) & 1
            if b12 == 0 and rs2 == 0: return "jump", 2
            if b12 == 1 and rs2 == 0 and rs1 != 0: return "jump", 2
        return "alu", 2
    w = h | (img[off + 2] << 16) | (img[off + 3] << 24); op = w & 0x7F
    if op == 0x03: return "load", 4
    if op == 0x23: return "store", 4
    if op == 0x63: return "branch", 4
    if op in (0x6F, 0x67): return "jump", 4
    if op == 0x33 and (w >> 25) == 1: return ("mul" if ((w >> 12) & 7) < 4 else "div"), 4
    return "alu", 4


cls_cache = {}
cnt = {k: 0 for k in ("alu", "load", "store", "branch", "branch_taken", "jump", "mul", "div")}
st = {"prev": None}


def hook(uc_, addr, size, _):
    c = cls_cache.get(addr)
    if c is None: c = classify(addr)[0]; cls_cache[addr] = c
    p = st["prev"]
    if p is not None and p[0] == "branch" and addr != p[1]: cnt["branch_taken"] += 1
    cnt[c] += 1; st["prev"] = (c, addr + size)


uc.hook_add(UC_HOOK_CODE, hook, begin=FLASH, end=FLASH + FLASH_SZ - 0x2000)


def call(name, *args):
    for k in cnt: cnt[k] = 0
    st["prev"] = None
    for r, a in zip((UC_RISCV_REG_A0, UC_RISCV_REG_A1, UC_RISCV_REG_A2), args): uc.reg_write(r, a)
    uc.reg_write(UC_RISCV_REG_SP, RAM + RAM_SZ); uc.reg_write(UC_RISCV_REG_RA, HALT)
    uc.emu_start(sym[name], HALT)
    snap = dict(cnt)
    return uc.reg_read(UC_RISCV_REG_A0), snap


tv = read_tv(f"fedelzeti/tv{N}.txt")
sb = bytearray(M1 * NMAX)
for i in range(M1):
    for t in range(tv["n"]): sb[i * NMAX + t] = tv["sign"][i][t] & 0xFF
uc.mem_write(sym["SIGN"], bytes(sb))

total = lambda c: c["alu"] + c["load"] + c["store"] + c["branch"] + c["jump"] + c["mul"] + c["div"]
# feltevesek (ciklus/utasitas; ESP32-C3 4 lepcsos in-order pipeline; a pontos ertekek NEM ismertek -> tartomany)
OPT = dict(alu=1, load=1, store=1, branch=1, branch_taken=0, jump=1, mul=1, div=20)       # optimista
PES = dict(alu=1, load=2, store=1, branch=1, branch_taken=2, jump=2, mul=4, div=35)       # pesszimista
cyc = lambda c, m: sum(c[k] * m[k] for k in m)

print(f"E27b | ESP32-C3 (RV32IMC, FPU nélkül, 160 MHz) emulátor | N = {tv['n']}")
print(f"kód: .text {code_size} B (+ .rodata {ro_size} B) | statikus RAM: .bss {bss_size} B + .data {data_size} B (benne a tesztkeret másolata: SIGN {M1*NMAX} B + a szűrő saját sign-táblája {M1*NMAX} B)")
call("c3_init", tv["n"])
init_total = 0
for s in tv["samples"]:
    _, c = call("c3_push", f2i(s[0]), f2i(s[1]), f2i(s[2])); init_total += total(c)
_, cd = call("c3_derive")
call("c3_dump")
out = struct.unpack("<8f", bytes(uc.mem_read(sym["OUT"], 32)))
print(f"\nELLENŐRZÉS a Python (f64) referenciával: theta emul {out[0]:.6f}, {out[1]:.6f} | Python {tv['th'][0]:.6f}, {tv['th'][1]:.6f} | max eltérés {max(abs(out[0]-tv['th'][0]), abs(out[1]-tv['th'][1])):.2e}")
print(f"  W emul {out[2]:.5f}, {out[3]:.5f}, {out[4]:.5f} | Python {tv['w'][0]:.5f}, {tv['w'][1]:.5f}, {tv['w'][2]:.5f} | nem-pd ellipszoid: emul {int(out[5])}, Python {tv['pd']}")
mx = 0.0
for h in tv["hx"]:
    r, _ = call("c3_hx", f2i(h[0]), f2i(h[1])); mx = max(mx, abs(i2f(r) - h[2]))
print(f"  burok-tartószám (16 irány): max abs eltérés {mx:.2e}")
_, cp = call("c3_push", f2i(tv["samples"][0][0]), f2i(tv["samples"][0][1]), f2i(tv["samples"][0][2]))
_, cd = call("c3_derive")
gates = []
dmax = 0.0; cnt_cmp = 0
for g in tv["gate"]:
    r, c = call("c3_gate", f2i(g[0]), f2i(g[1])); u = i2f(r); stt = int.from_bytes(uc.mem_read(sym["STATE"], 4), "little")
    gates.append((total(c), stt, c))
    if not math.isnan(g[2]): dmax = max(dmax, abs(u - g[2])); cnt_cmp += 1
print(f"  kapu (u): {cnt_cmp} eset, max eltérés {dmax:.4f} (a Python rács 0.025 lépésközű)")

gs = sorted(gates, key=lambda x: x[0])
names = {0: "átengedett", 1: "korlátozott", 2: "vész", 3: "tartalék"}
print("\nUTASÍTÁSSZÁM / LÉPÉS (MÉRT, emulátor):")
print(f"  push (1 minta, inkrementális): {total(cp)}   derive (θ, W, 99 ellipszoid): {total(cd)}")
print(f"  kapu: legkisebb {gs[0][0]} ({names[gs[0][1]]}), medián {gs[len(gs)//2][0]} ({names[gs[len(gs)//2][1]]}), LEGROSSZABB {gs[-1][0]} ({names[gs[-1][1]]})")
print(f"  (a kezdeti ablak-feltöltés N = {tv['n']} minta: {init_total} utasítás; ez csak felszállás előtt, egyszer)")
typ = (total(cp) + total(cd) + gs[len(gs) // 2][0]); wc = (total(cp) + total(cd) + gs[-1][0])
print(f"  EGY VEZÉRLÉSI LÉPÉS: tipikus {typ} utasítás, legrosszabb {wc} utasítás")

def mix(parts):
    return {k: sum(p[k] for p in parts) for k in parts[0]}
c_typ = mix([cp, cd, gs[len(gs) // 2][2]]); c_wc = mix([cp, cd, gs[-1][2]])
print("\nBECSÜLT IDŐ az ESP32-C3-on (utasításszám MÉRT; ciklus/utasítás FELTEVÉS: optimista = 1 ciklus, osztás 20; pesszimista = betöltés 2, vett ugrás +2, szorzás 4, osztás 35):")
for nm, c in (("tipikus", c_typ), ("legrosszabb", c_wc)):
    lo, hi = cyc(c, OPT), cyc(c, PES)
    print(f"  {nm:12s}: {lo/160e6*1e3:6.3f} – {hi/160e6*1e3:6.3f} ms  ({lo} – {hi} ciklus 160 MHz-en)")
print(f"  utasítás-összetétel (legrosszabb lépés): " + ", ".join(f"{k} {v}" for k, v in c_wc.items()))
lo, hi = cyc(c_wc, OPT) / 160e6, cyc(c_wc, PES) / 160e6
print(f"  => a legrosszabb lépés {hi*1e3:.2f} ms (pesszimista) -> {1/(hi*4):.0f} Hz-ig tartható 4x-es tartalékkal; 20 Hz-es (50 ms) hurok CPU-terhelése {hi/0.05*100:.1f}%")
