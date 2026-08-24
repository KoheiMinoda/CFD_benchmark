import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================
# CONFIG
# ============================================================
CSV_FILE   = "./force_trackingforce_tracking.csv"
OUTPUT_DIR = "./FORCE_FIGS"
os.makedirs(OUTPUT_DIR, exist_ok=True)

H  = 0.04
U0 = 0.535
T_shed = H / (0.132 * U0)   # 放出周期の目安 ≈ 0.566 s

T_CUT = 10.0   # この時刻より前を計算・描画から除く [s]

# ============================================================
# 1. CSV 読み込み → 最初の T_CUT 秒を除去
# ============================================================
data = np.loadtxt(CSV_FILE, delimiter=",", skiprows=1)
# 空白/タブ区切りなら: data = np.loadtxt(CSV_FILE, skiprows=1)

t_all = data[:, 0]
keep  = t_all >= T_CUT
data  = data[keep]

t   = data[:, 0]
Cd  = data[:, 4]
Cl  = data[:, 5]

print(f"after cut (t >= {T_CUT}s): {len(t)} steps, "
      f"t = {t[0]:.3f} ... {t[-1]:.3f} s "
      f"(~{(t[-1]-t[0])/T_shed:.1f} shedding periods)")

# ============================================================
# 2. 統計（全て 10秒以降のデータで計算される）
# ============================================================
Cd_mean = Cd.mean()
Cl_mean = Cl.mean()
Cl_rms  = np.sqrt(np.mean((Cl - Cl_mean)**2))

print(f"mean Cd = {Cd_mean:.3f}")
print(f"mean Cl = {Cl_mean:.3f},  Cl' (rms) = {Cl_rms:.3f}")

# ============================================================
# 3. St を Cl の FFT から算出
# ============================================================
cl = Cl - Cl_mean
dt = np.mean(np.diff(t))

fft   = np.abs(np.fft.rfft(cl * np.hanning(len(cl))))
freqs = np.fft.rfftfreq(len(cl), d=dt)
f_peak = freqs[np.argmax(fft[1:]) + 1]
St = f_peak * H / U0
print(f"St = {St:.3f}  (f_peak = {f_peak:.3f} Hz)")

# ============================================================
# 4. 時系列プロット（Cd, Cl）
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

ax1.plot(t, Cd, color="tab:blue", lw=0.8)
ax1.axhline(Cd_mean, color="k", ls="--", lw=0.8,
            label=rf"$\overline{{C_D}}$ = {Cd_mean:.3f}")
ax1.set_ylabel(r"$C_D$", fontsize=12)
ax1.legend(loc="upper right", fontsize=10)
ax1.grid(alpha=0.3)

ax2.plot(t, Cl, color="tab:red", lw=0.8)
ax2.axhline(Cl_mean, color="k", ls="--", lw=0.8)
ax2.set_ylabel(r"$C_L$", fontsize=12)
ax2.set_xlabel(r"$t$ [s]", fontsize=12)
ax2.legend([rf"$C_L'$ (rms) = {Cl_rms:.3f},  $St$ = {St:.3f}"],
           loc="upper right", fontsize=10)
ax2.grid(alpha=0.3)

fig.suptitle(rf"Force coefficients time history  ($t \geq$ {T_CUT} s)",
             fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "force_timeseries.png"),
            dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: force_timeseries.png")

# ============================================================
# 5. Cl のスペクトル
# ============================================================
fig, ax = plt.subplots(figsize=(7, 4))
St_axis = freqs * H / U0
ax.plot(St_axis, fft, color="tab:red", lw=1.0)
ax.axvline(St, color="k", ls="--", lw=0.8, label=rf"$St$ = {St:.3f}")
ax.set_xlabel(r"$St = fH/U_0$", fontsize=12)
ax.set_ylabel(r"PSD of $C_L$", fontsize=12)
ax.set_xlim(0, 1.0)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "cl_spectrum.png"),
            dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: cl_spectrum.png")

print("Done.")
