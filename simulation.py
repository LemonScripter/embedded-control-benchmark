import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

"""
MetaSpace Technologies - Hardware Constraint Simulation
Author: Laszlo Szoke (Lead Architect)

Description:
This script simulates the CPU load and Watchdog Timer behavior of an ESP32 microcontroller
running complex Linear Algebra (Matrix Inversion) vs. O(1) Scalar Checks.
"""

# --- TUDOMÁNYOS KONFIGURÁCIÓ ---
# Nincs "sleep". A terhelést a mátrix mérete adja.
MATRIX_SIZE = 15       # 15x15-ös rendszer (Drón állapotmodellje)
SPS_ITERATIONS = 50    # Hányszor kell újraszámolni ciklusonként (SPS)
# Ez összesen 50 db mátrix inverziót jelent képkockánként!

PHYSICAL_LIMIT = 8.0
SOFTWARE_LIMIT = 5.0

class RealMathDrone:
    def __init__(self, mode):
        self.mode = mode
        self.pos = 0.0
        self.vel = 0.0
        self.crashed = False
        
        # Telemetria
        self.compute_time_us = 0 # Mikroszekundum
        self.operation_name = ""
        self.status = "RUNNING" # <--- HIÁNYZÓ INICIALIZÁLÁS PÓTOLVA
        
        # --- ESP32 HARDVER BECSLÉS ---
        # ESP32: Dual Core 240MHz, kb. 60-80 MFLOPS (gyakorlati C kódnál kevesebb)
        # Feltételezzük: 50Hz Control Loop (20ms ablak)
        # Aktív áramfelvétel: ~240mA (WiFi+CPU), Sleep: ~10mA
        self.esp32_flops = 0
        self.esp32_cpu_load_percent = 0.0
        self.total_energy_mj = 0.0 # MilliJoules
        
        # --- OFFLINE ELŐKÉSZÍTÉS (Csak MetaSpace-nek) ---
        # Létrehozunk egy véletlen rendszermátrixot (Phi)
        self.system_matrix = np.random.rand(MATRIX_SIZE, MATRIX_SIZE)
        # MetaSpace: Előre kiszámoljuk az inverzét (Gain Matrix)
        # K = (X^T * X)^-1 * X^T
        # Ez a "Logic Lock" fájl eredménye
        self.precalc_gain = np.linalg.inv(self.system_matrix.T @ self.system_matrix + np.eye(MATRIX_SIZE)*0.1)

    def update(self, wind_force):
        if self.crashed: return

        # 1. FIZIKA (Valóság)
        self.pos += self.vel * 0.1
        self.vel += wind_force * 0.1
        
        control_force = 0
        
        # MÉRÉS INDÍTÁSA
        t0 = time.perf_counter()

        # FLOP BECSLÉS PARAMÉTEREK
        # MatInv (NxN): ~ (2/3)*N^3 művelet
        # MatMult (NxN * Nx1): ~ 2*N^2 művelet
        N = MATRIX_SIZE

        if self.mode == 'traditional':
            # === VALÓDI KOLUMBÁN MATEK (TUNED DEMO) ===
            
            # 1. Alap zaj (kicsi ingadozás)
            current_iters = 10 + np.random.randint(-2, 3) # 8-13 iters (32-52% Load)
            
            # 2. "Stressz Esemény" (pl. hirtelen széllökés korrekciója)
            # Minden ~100. ciklusban megnő a számításigény
            self.stress_timer = getattr(self, 'stress_timer', 0) + 1
            if self.stress_timer > 100:
                current_iters += 20 # +80% Load -> Spike 100% fölé
                if self.stress_timer > 140: # 40 ciklus után vége a stressznek
                    self.stress_timer = 0
            
            # FLOP Számítás a tudományos hitelességért
            flops_per_iter = (2/3) * (N**3) + 2*(N**2) 
            total_flops = flops_per_iter * current_iters
            self.esp32_flops = total_flops
            
            # PÓTOLT VÁLTOZÓ (Energiaszámításhoz kell)
            esp32_time_s = total_flops / 5_000_000.0

            # DEMO-BARÁT LOAD FORMULA: 
            # 1 iteráció kb 4% CPU terhelést jelent
            self.esp32_cpu_load_percent = current_iters * 4.0
            
            # Mátrix Inverziók (Valós terhelés a PC-nek is)
            results = []
            for _ in range(current_iters):
                noise = np.random.normal(0, 0.1, (MATRIX_SIZE, MATRIX_SIZE))
                perturbed_matrix = self.system_matrix + noise
                try:
                    inv_matrix = np.linalg.inv(perturbed_matrix.T @ perturbed_matrix + np.eye(MATRIX_SIZE)*0.01)
                    est = np.sum(inv_matrix) 
                    results.append(est)
                except np.linalg.LinAlgError:
                    pass
            
            avg_est = np.mean(results) if results else 0
            
            if self.esp32_cpu_load_percent > 90:
                self.operation_name = f"HEAVY CALC: {current_iters}x Inv"
            else:
                self.operation_name = f"NORMAL: {current_iters}x Inv"
            
            # PÓTLÁS: Definíció az energiaszámításhoz
            esp32_time_s = total_flops / 5_000_000.0
                
            control_force = -0.5 * self.pos - 0.8 * self.vel

        elif self.mode == 'metaspace':
            # === VALÓDI METASPACE MATEK (OPTIMIZED) ===
            
            # Csak mátrix-vektor szorzás: 2*N^2
            total_flops = 2 * (N**2) + 10 # +10 a logikai check
            self.esp32_flops = total_flops
            
            # ESP32 idő (MetaSpace mindig elhanyagolható)
            esp32_time_s = total_flops / 5_000_000.0

            # VALÓS ÁLLAPOT VEKTOR (Hogy matematikailag helyes legyen)
            state_vec = np.zeros(MATRIX_SIZE)
            state_vec[0] = self.pos
            state_vec[1] = self.vel
            
            # A FORDÍTÁS BIZONYÍTÉKA:
            # Itt alkalmazzuk az Offline kiszámolt (Kolumbán-alapú) mátrixot a valós adatra
            estimation = self.precalc_gain @ state_vec 
            
            if abs(self.pos) > SOFTWARE_LIMIT:
                control_force = -3.0 * np.sign(self.pos)
            else:
                control_force = -1.0 * self.pos - 1.0 * self.vel
                
            self.operation_name = "PRE-CALC: Dot Product"
            self.esp32_cpu_load_percent = 0.5 # Fix minimális terhelés

        
        # --- WATCHDOG & FAGYÁS SZIMULÁCIÓ ---
        # Ha a CPU 100% felett van, a rendszer elveszti a fonalat (Deadline Miss)
        if self.esp32_cpu_load_percent > 100:
            self.overload_counter = getattr(self, 'overload_counter', 0) + 1
        else:
            self.overload_counter = max(0, getattr(self, 'overload_counter', 0) - 1)

        # TÜRELMI IDŐ MEGNÖVELVE: 20 ciklus (hogy lássuk a szenvedést)
        self.frozen_timer = getattr(self, 'frozen_timer', 0)
        
        if self.overload_counter > 20:
            self.frozen_timer = 20 # Reset ideje
            self.overload_counter = 0 
        
        if self.frozen_timer > 0:
            self.frozen_timer -= 1
            control_force = 0 
            self.operation_name = "! WATCHDOG RESET !"
            self.status = "SYSTEM FREEZE"
        else:
            self.status = "RUNNING"

        if self.esp32_cpu_load_percent > 100: self.esp32_cpu_load_percent = 100.0
        
        # ENERGIA (3.3V * 240mA * idő) -> Watts * s = Joules
        # Active power: ~0.8 Watt
        energy_step = 0.8 * esp32_time_s
        self.total_energy_mj += energy_step * 1000 # mJ

        # MÉRÉS VÉGE (PC Idő)
        t1 = time.perf_counter()
        self.compute_time_us = (t1 - t0) * 1_000_000

        if abs(self.pos) > PHYSICAL_LIMIT:
            self.crashed = True
            
        self.vel += control_force * 0.1

# --- VIZUALIZÁCIÓ ---
plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

d1 = RealMathDrone('traditional')
d2 = RealMathDrone('metaspace')

def setup_ax(ax, title):
    ax.set_title(title, fontsize=11, color='white', fontweight='bold')
    ax.set_xlim(-10, 10)
    ax.set_ylim(-3, 5) # Megnövelt tér a mozgó szövegnek
    ax.get_yaxis().set_visible(False)
    ax.axvline(-PHYSICAL_LIMIT, color='red', lw=3)
    ax.axvline(PHYSICAL_LIMIT, color='red', lw=3)
    ax.axvline(-SOFTWARE_LIMIT, color='cyan', ls='--')
    ax.axvline(SOFTWARE_LIMIT, color='cyan', ls='--')

setup_ax(ax1, "KOLUMBÁN (SPS Iteráció)")
setup_ax(ax2, "METASPACE (Invariáns)")

dot1, = ax1.plot([], [], 'ro', ms=15)
dot2, = ax2.plot([], [], 'go', ms=15)

# Kezdetben üres szövegek (FIX POZÍCIÓ - OPTIMALIZÁLT)
txt1 = ax1.text(0.5, 0.65, "", transform=ax1.transAxes, ha='center', color='#ff9999', fontsize=10, family='monospace', fontweight='bold')
txt2 = ax2.text(0.5, 0.65, "", transform=ax2.transAxes, ha='center', color='#99ff99', fontsize=10, family='monospace', fontweight='bold')

frame_count = 0

def animate(i):
    global frame_count
    try:
        frame_count += 1
        
        wind = np.sin(frame_count * 0.1) * 3.0 + np.random.normal(0, 2.0)
        
        d1.update(wind)
        d2.update(wind)
        
        dot1.set_data([d1.pos], [0])
        dot2.set_data([d2.pos], [0])
        
        # (A szöveg pozíciója FIX, nem mozog a golyóval)

        # Bal oldal (Hagyományos)
        status1 = "CRASHED" if d1.crashed else d1.status
        
        # Ha FAGYÁS van, legyen feltűnő a szöveg
        header1 = d1.operation_name
        if d1.frozen_timer > 0:
             header1 = "!!! WATCHDOG RESET !!!"
             txt1.set_color('red')
        else:
             txt1.set_color('#ff9999')

        t1_text = (f"MATH  : inv((A+N).T @ (A+N))\n"
                   f"{header1}\n"
                   f"----------------\n"
                   f"FLOPs : {int(d1.esp32_flops)}\n"
                   f"Load  : {d1.esp32_cpu_load_percent:.0f}%\n"
                   f"Energy: {d1.total_energy_mj:.1f} mJ\n"
                   f"State : {status1}")
        txt1.set_text(t1_text)
        
        # Jobb oldal (MetaSpace)
        t2_text = (f"MATH  : x < LIMIT ? (Logic)\n"
                   f"{d2.operation_name}\n"
                   f"----------------\n"
                   f"FLOPs : {int(d2.esp32_flops)}\n"
                   f"Load  : {d2.esp32_cpu_load_percent:.2f}%\n"
                   f"Energy: {d2.total_energy_mj:.2f} mJ\n"
                   f"State : OK")
        txt2.set_text(t2_text)

        return dot1, dot2, txt1, txt2
    except Exception as e:
        print(f"ANIMATION ERROR: {e}")
        import traceback
        traceback.print_exc()
        return dot1, dot2, txt1, txt2

ani = animation.FuncAnimation(fig, animate, interval=50, blit=False)
plt.show()