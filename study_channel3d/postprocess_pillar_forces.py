import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===================== settings =====================

# attention to <time>
CSV_FILE = "./case_channal3d/RESU/20251117-1525/pillar_forces.csv"

# Use the last X% of samples for statistics (to avoid initial transients)
USE_LAST_FRACTION = 0.5 # e.g. 0.5 = last 50%, 0.2 = last 20%

# Geometry & reference values
D = 0.1 # characteristic length (pillar height or diameter) [m]
U_ref = 0.2 # reference velocity [m/s] (here: Re = 100 -> U_ref = 1.0)
nu = 1.0e-3 # kinematic viscosity [m^2/s] (rho=1, mu=1e-3)

# Reference values from the paper (experimental results)
REF_Cd_mean = 7.6
REF_Cl_mean = 0.07
REF_St = None
REF_dp = 0.17

# ========================================================

def main():
    df = pd.read_csv(CSV_FILE, skipinitialspace=True)

    # Expect columns: t, Fx, Fy, Fz, Cd, Cl, dp
    required = ["t", "Fx", "Fy", "Fz", "Cd", "Cl", "dp"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in {CSV_FILE}. "
                             f"Found columns: {list(df.columns)}")

    t_all = df["t"].to_numpy()
    Cd_all = df["Cd"].to_numpy()
    Cl_all = df["Cl"].to_numpy()
    dp_all = df["dp"].to_numpy()

    n_all = len(df)
    if n_all < 10:
        raise ValueError("Not enough samples in CSV (<10).")

    # ---------- select "statistical" window ----------
    i0 = int(n_all * (1.0 - USE_LAST_FRACTION))
    if i0 < 0:
        i0 = 0

    t = t_all[i0:]
    Cd = Cd_all[i0:]
    Cl = Cl_all[i0:]
    dp = dp_all[i0:]

    print("--------------------------------------------------")
    print(f"Total samples        : {n_all}")
    print(f"Using last           : {USE_LAST_FRACTION*100:.1f}%")
    print(f"Samples in window    : {len(t)}")
    print(f"Time window          : t = {t[0]:.6g}  ~  {t[-1]:.6g}")
    print("--------------------------------------------------")

    # ---------- basic statistics ----------
    Cd_mean = np.mean(Cd)
    Cd_std  = np.std(Cd)

    Cl_mean = np.mean(Cl)
    Cl_rms  = np.sqrt(np.mean((Cl - Cl_mean)**2))

    dp_mean = np.mean(dp)
    dp_std  = np.std(dp)

    # Re from U_ref, D, nu (for sanity check)
    Re = U_ref * D / nu

    print(f"Re (from U_ref, D, nu): {Re:.3f}")
    print("---- Time-averaged quantities (over selected window) ----")
    print(f"<Cd>   = {Cd_mean:.6f}  (std = {Cd_std:.6f})")
    print(f"<Cl>   = {Cl_mean:.6f}  (RMS around mean = {Cl_rms:.6f})")
    print(f"<dp>   = {dp_mean:.6f}  (std = {dp_std:.6f})")

    # ---------- Strouhal number from Cl(t) ----------
    # Assume almost constant dt
    dt_array = np.diff(t)
    dt = np.mean(dt_array)
    if dt <= 0.0:
        raise ValueError("Non-positive dt detected.")

    fs = 1.0 / dt     # sampling frequency
    N = len(Cl)

    # Remove mean before FFT
    Cl_0mean = Cl - Cl_mean

    # real FFT
    freqs = np.fft.rfftfreq(N, d=dt)
    fft_vals = np.fft.rfft(Cl_0mean)
    amp = np.abs(fft_vals)

    # Ignore zero frequency
    # and look for peak in a reasonable band (e.g. f > 0)
    # Optionally we can restrict to f < fs/2
    mask = freqs > 0.0
    freqs_pos = freqs[mask]
    amp_pos = amp[mask]

    if len(freqs_pos) > 0:
        j_peak = np.argmax(amp_pos)
        f_peak = freqs_pos[j_peak]
        St = f_peak * D / U_ref

        print("---- Spectral analysis of Cl(t) ----")
        print(f"dt       = {dt:.6e}")
        print(f"fs       = {fs:.6e}")
        print(f"N (FFT)  = {N}")
        print(f"f_peak   = {f_peak:.6e}  [Hz]")
        print(f"St       = {St:.6f}")
    else:
        f_peak = np.nan
        St = np.nan
        print("No positive frequencies found for FFT (too few samples?).")

    # ---------- Compare with reference values (if provided) ----------
    def rel_err(num, ref):
        return 100.0 * (num - ref) / ref if ref is not None and ref != 0.0 else None

    print("--------------------------------------------------")
    if REF_Cd_mean is not None:
        err = rel_err(Cd_mean, REF_Cd_mean)
        print(f"Ref <Cd> = {REF_Cd_mean:.6f},  error = {err:.2f} %")
    else:
        print("Ref <Cd> : (not set; fill REF_Cd_mean in script if you wish)")

    if REF_Cl_mean is not None:
        err = rel_err(Cl_mean, REF_Cl_mean)
        print(f"Ref <Cl> = {REF_Cl_mean:.6f},  error = {err:.2f} %")
    else:
        print("Ref <Cl> : (not set; fill REF_Cl_mean in script if you wish)")

    if REF_St is not None and not np.isnan(St):
        err = rel_err(St, REF_St)
        print(f"Ref St   = {REF_St:.6f},      error = {err:.2f} %")
    else:
        print("Ref St   : (not set; fill REF_St in script if you wish)")

    if REF_dp is not None:
        err = rel_err(dp_mean, REF_dp)
        print(f"Ref <dp> = {REF_dp:.6f},  error = {err:.2f} %")
    else:
        print("Ref <dp> : (not set; fill REF_dp in script if you wish)")
    print("--------------------------------------------------")

    # ---------- Plots ----------
    fig1, ax1 = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    ax1[0].plot(t, Cd, label="Cd")
    ax1[0].set_ylabel("Cd")
    ax1[0].grid(True)

    ax1[1].plot(t, Cl, label="Cl")
    ax1[1].set_ylabel("Cl")
    ax1[1].grid(True)

    ax1[2].plot(t, dp, label="dp")
    ax1[2].set_ylabel("dp")
    ax1[2].set_xlabel("time")
    ax1[2].grid(True)

    fig1.suptitle("Time series (window used for statistics)")
    fig1.tight_layout()

    # Spectrum of Cl
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 4))
    ax2.plot(freqs_pos, amp_pos, "-")
    ax2.set_xlim(0, max(freqs_pos) if len(freqs_pos) > 0 else 1.0)
    ax2.set_xlabel("frequency [Hz]")
    ax2.set_ylabel("|Cl_hat(f)|")
    ax2.set_title("Cl spectrum")
    ax2.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
