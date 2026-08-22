import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from scipy.interpolate import griddata

# ============================================================
# CONFIG
# ============================================================
CASE_FILE  = "../POSTS/RESULTS_FLUID_DOMAIN.case"
OUTPUT_DIR = "./POD_FIGS"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SNAP_START = 200
SNAP_END   = 600 # None にすると最終スナップまで
SNAP_STRIDE = 1

SLICE_NORMAL = (0, 0, 1)
SLICE_ORIGIN = (0, 0, 0.002)

N_MODES = 6
GRID_NX = 600
GRID_NY = 400

H = 0.04
U0 = 0.535 # 基準速度（無次元化用）
CYL_REGION = (-0.5, 0.5, -0.5, 0.5) # x/H 単位で柱を描く

# 描画範囲（x/H, y/H 単位）
XLIM_H = (-2, 15)
YLIM_H = (-4, 4)

# ============================================================
# HELPER
# ============================================================
def draw_cylinder(ax):
    x0, x1, y0, y1 = CYL_REGION
    ax.add_patch(plt.Rectangle((x0, y0), x1-x0, y1-y0,
                 linewidth=0.8, edgecolor="black",
                 facecolor="lightgray", zorder=5))

# ============================================================
# スナップショット読み込み（u, v 成分）
# ============================================================
reader = pv.get_reader(CASE_FILE)
time_values = np.array(reader.time_values)

end = SNAP_END if SNAP_END is not None else len(time_values)
snap_indices = list(range(SNAP_START, end, SNAP_STRIDE))
M = len(snap_indices)

# サンプリング間隔 時間係数のFFT用
dt_snap = np.mean(np.diff(time_values[snap_indices]))
print(f"{M} snapshots, dt_snap = {dt_snap:.4e} s")

def get_sliced_uv(idx):
    t = float(time_values[idx])
    reader.set_active_time_value(t)
    mesh = reader.read()
    fluid = mesh["Fluid domain"] if isinstance(mesh, pv.MultiBlock) else mesh
    fluid_pt = fluid.cell_data_to_point_data()
    vel = fluid_pt.point_data["Velocity"]
    fluid_pt["u_vel"] = vel[:, 0]
    fluid_pt["v_vel"] = vel[:, 1]
    sliced = fluid_pt.slice(normal=SLICE_NORMAL, origin=SLICE_ORIGIN)
    pts = np.array(sliced.points)
    return pts, np.array(sliced["u_vel"]), np.array(sliced["v_vel"])

# 共通グリッド x/H, y/H 単位で作る
xi = np.linspace(XLIM_H[0]*H, XLIM_H[1]*H, GRID_NX)
yi = np.linspace(YLIM_H[0]*H, YLIM_H[1]*H, GRID_NY)
Xi, Yi = np.meshgrid(xi, yi)
Xi_H, Yi_H = Xi/H, Yi/H # 描画用の無次元座標

def interp(pts, vals):
    return griddata((pts[:,0], pts[:,1]), vals, (Xi, Yi), method="linear")

Ngrid = GRID_NX * GRID_NY
# u, v を縦に積む: X は (2*Ngrid, M)
X = np.zeros((2*Ngrid, M))
for j, idx in enumerate(snap_indices):
    pts, uu, vv = get_sliced_uv(idx)
    X[:Ngrid, j]      = interp(pts, uu).ravel()
    X[Ngrid:2*Ngrid, j] = interp(pts, vv).ravel()
    if (j+1) % 20 == 0:
        print(f"  loaded {j+1}/{M}")

X = np.nan_to_num(X, nan=0.0)

# ============================================================
# 平均場と変動場
# ============================================================
X_mean = X.mean(axis=1, keepdims=True)
Xp = X - X_mean

u_mean = X_mean[:Ngrid, 0].reshape(GRID_NY, GRID_NX)
v_mean = X_mean[Ngrid:2*Ngrid, 0].reshape(GRID_NY, GRID_NX)

# ============================================================
# Snapshot POD
# ============================================================
C = Xp.T @ Xp
lam, Psi = np.linalg.eigh(C)
order = np.argsort(lam)[::-1]
lam = lam[order]; Psi = Psi[:, order]
lam[lam < 0] = 0.0
sigma = np.sqrt(lam)

Phi = np.zeros((2*Ngrid, M))
for k in range(M):
    if sigma[k] > 1e-12:
        Phi[:, k] = (Xp @ Psi[:, k]) / sigma[k]

energy = lam / np.sum(lam)

# 時間係数 a_k(t) = sigma_k * psi_k
A = (np.diag(sigma) @ Psi.T)     # (M, M): 行 k が a_k(t)

# 各モードの St を時間係数のFFTで算出
def mode_St(ak):
    ak = ak - ak.mean()
    fft = np.abs(np.fft.rfft(ak * np.hanning(len(ak))))
    freqs = np.fft.rfftfreq(len(ak), d=dt_snap)
    f_peak = freqs[np.argmax(fft[1:]) + 1]   # DC除く
    return f_peak * H / U0

St_modes = [mode_St(A[k]) for k in range(N_MODES)]
print("Mode energy (%):", np.round(energy[:N_MODES]*100, 2))
print("Mode St:", np.round(St_modes, 3))

# ============================================================
# 描画：平均場 + 各モードの φx, φy を2列で
# ============================================================
fig, axes = plt.subplots(N_MODES+1, 2, figsize=(14, 3.2*(N_MODES+1)))
fig.suptitle(r"POD modes of $U$  (normalised: $\langle\phi_i,\phi_i\rangle=1$)",
             fontsize=14, y=0.995)

def plot_field(ax, field, title, clim, cmap="RdBu_r", center0=True):
    if center0:
        norm = mcolors.TwoSlopeNorm(vmin=-clim, vcenter=0.0, vmax=clim)
        levels = np.linspace(-clim, clim, 31)
    else:
        norm = None
        levels = np.linspace(clim[0], clim[1], 31)
    cf = ax.contourf(Xi_H, Yi_H, field, levels=levels, cmap=cmap,
                     norm=norm, extend="both")
    draw_cylinder(ax)
    ax.set_aspect("equal")
    ax.set_xlim(XLIM_H); ax.set_ylim(YLIM_H)
    ax.set_ylabel(r"$y/H$")
    ax.set_title(title, fontsize=10)
    plt.colorbar(cf, ax=ax, fraction=0.025, pad=0.01)

# 平均場 (u/U0, v/U0)
plot_field(axes[0,0], u_mean/U0, r"time mean  $\bar{u}/U_0$",
           clim=(-0.2, 1.6), cmap="jet", center0=False)
plot_field(axes[0,1], v_mean/U0, r"time mean  $\bar{v}/U_0$", clim=0.8)

# 各モード
for k in range(N_MODES):
    phi_u = Phi[:Ngrid, k].reshape(GRID_NY, GRID_NX)
    phi_v = Phi[Ngrid:2*Ngrid, k].reshape(GRID_NY, GRID_NX)
    
    amax = max(np.nanmax(np.abs(phi_u)), np.nanmax(np.abs(phi_v)))
    scale = amax if amax > 0 else 1.0
    ttl = (rf"mode {k+1}  $\phi_x$   "
           rf"$\lambda_i/\Sigma\lambda$ = {energy[k]*100:.2f}%,  "
           rf"$St$ = {St_modes[k]:.3f}")
    ttl_v = ttl.replace(r"$\phi_x$", r"$\phi_y$")
    clim = 4.0   # 画像に合わせた固定スケール（要調整）
    plot_field(axes[k+1,0], phi_u/scale*clim, ttl, clim=clim)
    plot_field(axes[k+1,1], phi_v/scale*clim, ttl_v, clim=clim)

axes[-1,0].set_xlabel(r"$x/H$")
axes[-1,1].set_xlabel(r"$x/H$")

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "pod_U_modes.png"),
            dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: pod_U_modes.png")
