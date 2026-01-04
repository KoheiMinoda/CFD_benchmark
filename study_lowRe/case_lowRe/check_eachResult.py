import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# ==========================================
# setting
# ==========================================
csv_file_path = 'Re41/DragDecompositiondrag_monitoring.csv' 
output_image  = 'Re10_drag_history.png'
TRIM_START_RATIO = 0.01

# 解析用：後半何割を使って平均値（定常値）を算出するか
CALC_AVG_RATIO = 0.2 

# グラフの見た目設定
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['font.size'] = 12

# ==========================================
# Main Process
# ==========================================
def main():
    if not os.path.exists(csv_file_path):
        print(f"[ERROR]: Can not find '{csv_file_path}'")
        return

    try:
        # --- 1. Load Data ---
        df = pd.read_csv(csv_file_path)
        
        # Column cleaning
        df.columns = [c.strip().replace('"', '').replace('#', '').strip() for c in df.columns]

        # Identify columns
        time_col = next((c for c in df.columns if c.lower() in ['t', 'time', 'physical_time']), None)
        cd_p_col = next((c for c in df.columns if 'cd_pressure' in c.lower()), None)
        cd_v_col = next((c for c in df.columns if 'cd_viscous' in c.lower()), None)
        cd_t_col = next((c for c in df.columns if 'cd_total' in c.lower()), None)
        cl_t_col = next((c for c in df.columns if 'cl_total' in c.lower()), None)
        re_col   = next((c for c in df.columns if 're' in c.lower()), None)

        if not (time_col and cd_p_col and cd_v_col):
            print("[ERROR]: Essential columns not found.")
            return

        # --- 2. Extract Data ---
        t_data = df[time_col].values
        cd_p   = df[cd_p_col].values
        cd_v   = df[cd_v_col].values
        cd_t   = df[cd_t_col].values if cd_t_col else (cd_p + cd_v)
        cl_t   = df[cl_t_col].values if cl_t_col else np.zeros_like(t_data)
        
        # Get Reynolds number (from the end of data)
        current_Re = df[re_col].values[-1] if re_col else 0.0

        # ==========================================
        # 【変更点】Trim Start Data (最初の1%をカット)
        # ==========================================
        n_total = len(t_data)
        n_cut   = int(n_total * TRIM_START_RATIO)
        
        if n_cut > 0 and n_cut < n_total:
            print(f"Trimming start: cut first {n_cut} steps ({TRIM_START_RATIO*100}%)")
            t_data = t_data[n_cut:]
            cd_p   = cd_p[n_cut:]
            cd_v   = cd_v[n_cut:]
            cd_t   = cd_t[n_cut:]
            cl_t   = cl_t[n_cut:]
        
        # --- 3. Calculate Steady State Values ---
        # (現在残っているデータの) 最後の N% を使用
        n_points = len(t_data)
        n_start_avg  = int(n_points * (1.0 - CALC_AVG_RATIO))
        
        if n_start_avg < n_points:
            avg_cd_p = np.mean(cd_p[n_start_avg:])
            avg_cd_v = np.mean(cd_v[n_start_avg:])
            avg_cd_t = np.mean(cd_t[n_start_avg:])
            
            print("\n" + "="*40)
            print(f" ANALYSIS RESULT (Re = {current_Re:.2f})")
            print("="*40)
            print(f" Pressure Drag (Cd_p) : {avg_cd_p:.5f}")
            print(f" Viscous  Drag (Cd_v) : {avg_cd_v:.5f}")
            print(f" Total    Drag (Cd_t) : {avg_cd_t:.5f}")
            print("="*40 + "\n")

        # --- 4. Plotting ---
        fig, ax1 = plt.subplots(figsize=(10, 6))

        ax1.plot(t_data, cd_t, color='black',     linestyle='-',  linewidth=2.0, label=r'$C_{D,total}$')
        ax1.plot(t_data, cd_p, color='tab:red',   linestyle='--', linewidth=1.5, label=r'$C_{D,pressure}$')
        ax1.plot(t_data, cd_v, color='tab:blue',  linestyle='--', linewidth=1.5, label=r'$C_{D,viscous}$')
        
        # Lift check
        if np.max(np.abs(cl_t)) > 0.001:
             ax2 = ax1.twinx()
             ax2.plot(t_data, cl_t, color='tab:green', linestyle=':', linewidth=1.0, alpha=0.7, label=r'$C_{L}$')
             ax2.set_ylabel(r'Lift Coefficient $C_L$', color='tab:green', fontsize=14)
             ax2.tick_params(axis='y', labelcolor='tab:green')
             # Legend hack
             ax1.plot([], [], color='tab:green', linestyle=':', label=r'$C_{L}$')

        ax1.set_xlabel(r'Physical Time [s]', fontsize=14)
        ax1.set_ylabel(r'Drag Coefficient $C_D$', fontsize=14)
        ax1.set_title(f'Drag History at Re = {current_Re:.1f}', fontsize=16)
        
        ax1.grid(True, which='major', linestyle=':', color='gray', alpha=0.5)
        ax1.legend(loc='center right', fontsize=12, frameon=True, edgecolor='black', fancybox=False)

        plt.tight_layout()
        plt.savefig(output_image, dpi=300)
        print(f"[COMPLETE] Saved plot as '{output_image}'")
        plt.show()

    except Exception as e:
        print(f"[ERROR] unexpected error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
