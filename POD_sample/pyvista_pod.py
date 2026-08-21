import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import os
from scipy.interpolate import griddata

# ============================================================
# CONFIG
# ============================================================
CASE_FILE  = "../POSTS/RESULTS_FLUID_DOMAIN.case"   # POD用出力を指す
OUTPUT_DIR = "./POD_FIGS"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# POD に使うスナップショットの時刻インデックス範囲
# 初期過渡を除き、リミットサイクルに乗った区間を使う
SNAP_START = 100      # このステップ以降を使用
SNAP_END   = None     # None なら最後まで
SNAP_STRIDE = 1       # 間引き（1=全部）

SLICE_NORMAL = (0, 0, 1)
SLICE_ORIGIN = (0, 0, 0.0628)   # z = Lz/2 = pi*H/2 ≈ 0.0628（スパン中央で切る）

N_MODES = 6           # 可視化する上位モード数

GRID_NX = 600
GRID_NY = 400

# 柱の位置（Case043: 原点中心, 一辺 H=0.04）
H = 0.04
CYL_REGION = (-H/2, H/2, -H/2, H/2)   # (xmin, xmax, ymin, ymax)

# 描画範囲（後流を見る）
XLIM = (-2*H, 8*H)
YLIM = (-3*H, 3*H)

# ============================================================
# HELPER: 柱を塗る
# ============================================================
def draw_cylinder(ax, region, color="dimgray", zorder=5):
    x_min, x_max, y_min, y_max = region
    rect = plt.Rectangle(
        (x_min, y_min), x_max - x_min, y_max - y_min,
        linewidth=0.8, edgecolor="black", facecolor=color, zorder=zorder
    )
    ax.add_patch(rect)

# ============================================================
# 1. 全スナップショットをスライスして u 成分を行列化
# ============================================================
reader = pv.get_reader(CASE_FILE)
time_values = np.array(reader.time_values)

end = SNAP_END if SNAP_END is not None else len(time_values)
snap_indices = list(range(SNAP_START, end, SNAP_STRIDE))
M = len(snap_indices)
print(f"{M} snapshots, t = {time_values[snap_indices[0]]:.3f} ... "
      f"{time_values[snap_indices[-1]]:.3f} s")

# 最初のスナップショットで共通グリッドを定義
def get_sliced_u(idx):
    t = float(time_values[idx])
    reader.set_active_time_value(t)
    mesh = reader.read()
    fluid = mesh["Fluid domain"] if isinstance(mesh, pv.MultiBlock) else mesh
    fluid_pt = fluid.cell_data_to_point_data()   # セル→点データ（添付コードと同じ）
    u = fluid_pt.point_data["Velocity"][:, 0]    # u成分だけ
    fluid_pt["u_vel"] = u
    sliced = fluid_pt.slice(normal=SLICE_NORMAL, origin=SLICE_ORIGIN)
    pts = np.array(sliced.points)
    uvals = np.array(sliced["u_vel"])
    return pts, uvals

# 共通の補間グリッドを作る（最初のスライスの点群範囲で固定）
pts0, u0 = get_sliced_u(snap_indices[0])
x2d, y2d = pts0[:, 0], pts0[:, 1]
xi = np.linspace(XLIM[0], XLIM[1], GRID_NX)
yi = np.linspace(YLIM[0], YLIM[1], GRID_NY)
Xi, Yi = np.meshgrid(xi, yi)

def interp_to_grid(pts, vals):
    return griddata((pts[:,0], pts[:,1]), vals, (Xi, Yi), method="linear")

# スナップショット行列 X: (Ngrid, M)
Ngrid = GRID_NX * GRID_NY
X = np.zeros((Ngrid, M))
for j, idx in enumerate(snap_indices):
    pts, uvals = get_sliced_u(idx)
    Zi = interp_to_grid(pts, uvals)
    X[:, j] = Zi.ravel()
    if (j+1) % 20 == 0:
        print(f"  loaded {j+1}/{M}")

# NaN（グリッド外/柱内部）を0で埋める（補間の穴対策）
X = np.nan_to_num(X, nan=0.0)

# ============================================================
# 2. 平均を引いて変動場に
# ============================================================
X_mean = X.mean(axis=1, keepdims=True)
Xp = X - X_mean

# ============================================================
# 3. Snapshot POD（method of snapshots: M×M で軽い）
# ============================================================
# 相関行列 C = Xp^T Xp  (M×M)
C = Xp.T @ Xp
lam, Psi = np.linalg.eigh(C)          # 昇順固有値
idx_sort = np.argsort(lam)[::-1]      # 降順に
lam = lam[idx_sort]
Psi = Psi[:, idx_sort]
lam[lam < 0] = 0.0                     # 数値誤差の負値をクリップ
sigma = np.sqrt(lam)

# 空間モード Phi = Xp Psi / sigma
Phi = np.zeros((Ngrid, M))
for k in range(M):
    if sigma[k] > 1e-12:
        Phi[:, k] = (Xp @ Psi[:, k]) / sigma[k]

energy = lam / np.sum(lam)
print("Mode energy (%):", np.round(energy[:N_MODES]*100, 2))
print("Cumulative (%):", np.round(np.cumsum(energy[:N_MODES])*100, 2))

# ============================================================
# 4. 各モードの空間構造を PNG に（モード1..N_MODES）
# ============================================================
for k in range(N_MODES):
    phi_k = Phi[:, k].reshape(GRID_NY, GRID_NX)

    # モードの符号・振幅を正規化（見やすさのため最大絶対値で規格化）
    amax = np.nanmax(np.abs(phi_k))
    if amax > 0:
        phi_k = phi_k / amax

    fig, ax = plt.subplots(figsize=(8, 4))
    cf = ax.contourf(
        Xi, Yi, phi_k,
        levels=np.linspace(-1, 1, 21),
        cmap="RdBu_r",
        norm=mcolors.TwoSlopeNorm(vmin=-1, vcenter=0.0, vmax=1),
        extend="both", zorder=1
    )
    ax.contour(Xi, Yi, phi_k, levels=np.linspace(-1, 1, 11),
               colors="k", linewidths=0.3, alpha=0.4, zorder=2)

    draw_cylinder(ax, CYL_REGION, zorder=5)

    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label(r"$\phi_u$ (normalized)", fontsize=11)

    ax.set_xlabel("x [m]", fontsize=12)
    ax.set_ylabel("y [m]", fontsize=12)
    ax.set_title(
        rf"POD mode {k+1}  ($u$)  —  E = {energy[k]*100:.1f}%",
        fontsize=13
    )
    ax.set_aspect("equal")
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)

    fname = os.path.join(OUTPUT_DIR, f"pod_mode_{k+1:02d}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fname}")

# ============================================================
# 5. エネルギースペクトルも出す
# ============================================================
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(np.arange(1, N_MODES+1), energy[:N_MODES]*100, color="steelblue")
ax.set_xlabel("Mode number", fontsize=12)
ax.set_ylabel("Energy fraction [%]", fontsize=12)
ax.set_title("POD energy spectrum", fontsize=13)
for i in range(N_MODES):
    ax.text(i+1, energy[i]*100, f"{energy[i]*100:.1f}",
            ha="center", va="bottom", fontsize=9)
fig.savefig(os.path.join(OUTPUT_DIR, "pod_energy_spectrum.png"),
            dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: pod_energy_spectrum.png")

print("Done.")
