"""
Independent Python re-check of the two-story building models (no MATLAB needed).

What it does
  1. Builds the mechanical and electrical state-space models (same matrices as matlab/*.m).
  2. Re-solves both from their *second-order* equations (M x'' + C x' + K x = F, and the
     electrical mesh form L q'' + R q' + S q = V) with a different integrator.
  3. Prints peak values, the analogy residuals and the damping-sensitivity table.
  4. Writes results/metrics.md and the PNG figures in results/figures/.

Run:  python python/verify_and_plot.py      (needs numpy, scipy, matplotlib)
"""
from pathlib import Path
import numpy as np
from scipy import signal
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- parameters (G = 4)
G = 4
m1, m2 = 1000 * G, 800 / G            # kg      -> 4000, 200
k1, k2 = 1.0e6 * G, 0.8e6 / G         # N/m     -> 4e6, 2e5
c1, c2 = 1.0e6 * G, 8.0e6 / G         # N s/m   -> 4e6, 2e6
F0, w = 5000.0, 2 * np.pi * G         # N, rad/s  (4 Hz)
T_END, N = 5.0, 5001
t = np.linspace(0.0, T_END, N)
u = F0 * np.sin(w * t)


def mech_A(c1, c2):
    return np.array([[0, 1, 0, 0],
                     [-(k1 + k2) / m1, -(c1 + c2) / m1, k2 / m1, c2 / m1],
                     [0, 0, 0, 1],
                     [k2 / m2, c2 / m2, -k2 / m2, -c2 / m2]])


B = np.array([[0], [1 / m1], [0], [0]])

# force-voltage (Maxwell) analogy: m -> L, c -> R, k -> 1/C, F -> V, x -> q, v -> i
L1, L2, R1, R2, C1, C2 = m1, m2, c1, c2, 1 / k1, 1 / k2
A_e = np.array([[0, 1, 0, 0],
                [-(1 / (L1 * C1) + 1 / (L1 * C2)), -(R1 + R2) / L1, 1 / (L1 * C2), R2 / L1],
                [0, 0, 0, 1],
                [1 / (L2 * C2), R2 / L2, -1 / (L2 * C2), -R2 / L2]])


def lsim(A, Bm, tt=t, uu=u):
    s = signal.StateSpace(A, Bm, np.eye(4), np.zeros((4, 1)))
    return signal.lsim(s, uu, tt)[1]


y_m = lsim(mech_A(c1, c2), B)       # [x1, v1, x2, v2]
y_e = lsim(A_e, B)                  # [q1, i1, q2, i2]


# ---- independent second-order solves (different formulation + different integrator)
def solve_second_order(Mm, Cm, Km, force):
    Minv = np.linalg.inv(Mm)

    def f(tt, z):
        x, v = z[:2], z[2:]
        a = Minv @ (np.array([force(tt), 0.0]) - Cm @ v - Km @ x)
        return np.r_[v, a]

    sol = solve_ivp(f, (0, T_END), np.zeros(4), t_eval=t, method="Radau",
                    rtol=1e-10, atol=1e-16)
    return sol.y.T        # [x1, x2, v1, v2]


M = np.diag([m1, m2])
Cm = np.array([[c1 + c2, -c2], [-c2, c2]])
Km = np.array([[k1 + k2, -k2], [-k2, k2]])
so_m = solve_second_order(M, Cm, Km, lambda tt: F0 * np.sin(w * tt))
# electrical mesh equations: L q'' + R q' + (1/C) q = V, written independently from the C's
Le = np.diag([L1, L2])
Re = np.array([[R1 + R2, -R2], [-R2, R2]])
Se = np.array([[1 / C1 + 1 / C2, -1 / C2], [-1 / C2, 1 / C2]])
so_e = solve_second_order(Le, Re, Se, lambda tt: F0 * np.sin(w * tt))

x1, v1, x2, v2 = y_m.T
q1, i1, q2, i2 = y_e.T
drift_m, drift_e = np.abs(x1 - x2), np.abs(q1 - q2)


def rel(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b)) * 100


res = {
    "mech: state-space vs 2nd-order ODE   (x1, x2)": max(rel(x1, so_m[:, 0]), rel(x2, so_m[:, 1])),
    "mech: state-space vs 2nd-order ODE   (v1, v2)": max(rel(v1, so_m[:, 2]), rel(v2, so_m[:, 3])),
    "elec: state-space vs mesh-equation ODE (q1, q2)": max(rel(q1, so_e[:, 0]), rel(q2, so_e[:, 1])),
    "elec: state-space vs mesh-equation ODE (i1, i2)": max(rel(i1, so_e[:, 2]), rel(i2, so_e[:, 3])),
    "analogy: x vs q (mech SS vs elec SS)": max(rel(x1, q1), rel(x2, q2)),
    "analogy: v vs i (mech SS vs elec SS)": max(rel(v1, i1), rel(v2, i2)),
    "analogy: mech SS vs elec mesh-ODE (x1 vs q1)": rel(x1, so_e[:, 0]),
}
ss = t >= 3.0     # steady-state window (last 2 s)
peaks = {
    "peak x1 (abs)": np.max(np.abs(x1)), "peak x2 (abs)": np.max(np.abs(x2)),
    "peak v1 (abs)": np.max(np.abs(v1)), "peak v2 (abs)": np.max(np.abs(v2)),
    "peak inter-story drift": drift_m.max(), "steady-state drift (t>=3 s)": drift_m[ss].max(),
    "steady-state x1 amplitude (t>=3 s)": np.abs(x1[ss]).max(),
}

# ---------------------------------------------------------------- sensitivity
cases = {"Case 1\nc1 +25%": (1.25 * c1, c2), "Case 2\nc2 −25%": (c1, 0.75 * c2),
         "Case 3\nboth +25%": (1.25 * c1, 1.25 * c2), "Case 4\nboth −25%": (0.75 * c1, 0.75 * c2)}
base = drift_m.max()
sens = {}
sens_traces = {}
for name, (a, b) in cases.items():
    yy = lsim(mech_A(a, b), B)
    sens_traces[name] = yy
    d = np.abs(yy[:, 0] - yy[:, 2]).max()
    sens[name.replace("\n", " ")] = (d, (d - base) / base * 100)

eig = np.linalg.eigvals(mech_A(c1, c2))

# ---------------------------------------------------------------- console + markdown
print("Eigenvalues of A (rad/s):", np.round(eig.real, 4))
for k, v in peaks.items():
    print(f"{k:34s} {v:.4e}")
for k, v in res.items():
    print(f"{k:50s} {v:.3e} %")
for k, (d, p) in sens.items():
    print(f"{k:20s} max drift {d:.4e} m  change {p:+.2f} %")

md = ["# Computed metrics (generated by `python/verify_and_plot.py`)", "",
      "G = 4, F(t) = 5000 sin(8πt) N, 0–5 s, zero initial conditions.", "",
      "| Quantity | Value |", "|---|---:|"]
unit = {"x": "m", "v": "m/s", "d": "m"}
for k, v in peaks.items():
    md.append(f"| {k} | {v:.3e} {'m/s' if 'v1' in k or 'v2' in k else 'm'} |")
md += ["", "| Check | Max relative difference |", "|---|---:|"]
md += [f"| {k} | {v:.2e} % |" for k, v in res.items()]
md += ["", "| Damping case | Max inter-story drift (m) | Change vs. baseline |", "|---|---:|---:|",
       f"| Baseline | {base:.4e} | — |"]
md += [f"| {k} | {d:.4e} | {p:+.2f} % |" for k, (d, p) in sens.items()]
md += ["", "A eigenvalues (rad/s): " + ", ".join(f"{e:.4g}" for e in sorted(eig.real)), ""]
(ROOT / "results" / "metrics.md").write_text("\n".join(md), encoding="utf-8")

# ---------------------------------------------------------------- figures
BLUE, ORANGE = "#2a78d6", "#eb6834"          # validated categorical slots 1-2
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.titlecolor": INK, "axes.titleweight": "bold", "axes.titlesize": 11,
    "axes.titlelocation": "left", "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5,
    "legend.frameon": False, "lines.linewidth": 1.6,
})


def legend(ax, **kw):
    kw.pop("loc", None)
    leg = ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=kw.pop("ncol", 2), borderaxespad=0.2, **kw)
    for tx in leg.get_texts():
        tx.set_color(MUTED)


# Figure 1: overview ------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(12, 7.2), constrained_layout=True)
ax = axs[0, 0]
ax.plot(t, x1 * 1e6, color=BLUE, label="Floor 1  x₁")
ax.plot(t, x2 * 1e6, color=ORANGE, ls="--", label="Floor 2  x₂")
ax.set(title="Floor displacement", xlabel="Time (s)", ylabel="Displacement (µm)")
legend(ax, loc="upper right")
ax = axs[0, 1]
ax.plot(t, v1 * 1e3, color=BLUE, label="Floor 1  v₁")
ax.plot(t, v2 * 1e3, color=ORANGE, ls="--", label="Floor 2  v₂")
ax.set(title="Floor velocity", xlabel="Time (s)", ylabel="Velocity (mm/s)")
legend(ax, loc="upper right")
ax = axs[1, 0]
ax.plot(t, drift_m * 1e9, color=BLUE, lw=2.2, label="Mechanical")
ax.plot(t, drift_e * 1e9, color=ORANGE, ls=(0, (4, 3)), lw=1.6, label="Electrical")
ax.set(title="Drift, first 1.5 s", xlabel="Time (s)", ylabel="Drift (nm ≙ nC)", xlim=(0, 1.5))
legend(ax, loc="lower right")
ax = axs[1, 1]
names = list(sens.keys())
vals = [sens[n][1] for n in names]
bars = ax.bar(range(4), vals, color=BLUE, width=0.55)
ax.axhline(0, color=MUTED, lw=1)
ax.grid(axis="x", visible=False)
ax.set_xticks(range(4), ["c₁ +25%\n(case 1)", "c₂ −25%\n(case 2)", "c₁, c₂ +25%\n(case 3)", "c₁, c₂ −25%\n(case 4)"])
for b, vv in zip(bars, vals):
    ax.annotate(f"{vv:+.1f}%", (b.get_x() + b.get_width() / 2, vv), ha="center",
                va="bottom" if vv >= 0 else "top", xytext=(0, 4 if vv >= 0 else -4),
                textcoords="offset points", color=INK, fontsize=10, fontweight="bold")
ax.set(title="Damping sensitivity: change in peak drift", ylabel="Change vs. baseline (%)", ylim=(-55, 100))
fig.suptitle("Two-story building under a 4 Hz, 5 kN base excitation (G = 4)", x=0.01, ha="left",
             fontsize=14, fontweight="bold")
fig.savefig(FIG / "overview.png", dpi=170)
plt.close(fig)

# Figure 2: analogy + independent check ---------------------------------
fig, axs = plt.subplots(2, 2, figsize=(12, 6.2), constrained_layout=True)
pairs = [(x1 * 1e6, q1 * 1e6, "Floor 1: displacement x₁ ↔ charge q₁", "µm ≙ µC"),
         (x2 * 1e6, q2 * 1e6, "Floor 2: displacement x₂ ↔ charge q₂", "µm ≙ µC"),
         (v1 * 1e3, i1 * 1e3, "Floor 1: velocity v₁ ↔ current i₁", "mm/s ≙ mA"),
         (v2 * 1e3, i2 * 1e3, "Floor 2: velocity v₂ ↔ current i₂", "mm/s ≙ mA")]
for ax, (a, b, ttl, yl) in zip(axs.ravel(), pairs):
    ax.plot(t, a, color=BLUE, lw=2.2, label="Mechanical")
    ax.plot(t, b, color=ORANGE, ls=(0, (4, 3)), lw=1.4, label="Electrical analog")
    ax.set(title=ttl, xlabel="Time (s)", ylabel=yl)
    legend(ax, loc="upper right", ncol=2)
fig.suptitle("Force–voltage analogy: the circuit reproduces the building", x=0.01, ha="left",
             fontsize=14, fontweight="bold")
fig.savefig(FIG / "analogy_overlay.png", dpi=170)
plt.close(fig)

# Figure 3: sensitivity traces -------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(12, 6.4), constrained_layout=True, sharex=True, sharey=True)
lab = {"Case 1\nc1 +25%": "Case 1 · c₁ +25%", "Case 2\nc2 −25%": "Case 2 · c₂ −25%",
       "Case 3\nboth +25%": "Case 3 · c₁, c₂ +25%", "Case 4\nboth −25%": "Case 4 · c₁, c₂ −25%"}
for ax, (name, yy) in zip(axs.ravel(), sens_traces.items()):
    h0, = ax.plot(t, np.abs(y_m[:, 0] - y_m[:, 2]) * 1e9, color=MUTED, lw=1.1, ls="--", label="Baseline")
    h1, = ax.plot(t, np.abs(yy[:, 0] - yy[:, 2]) * 1e9, color=BLUE, lw=1.6, label="Modified damping")
    ax.set_title(lab[name])
    ax.text(0.99, 0.95, f"{sens[name.replace(chr(10), ' ')][1]:+.1f}% peak drift", transform=ax.transAxes,
            ha="right", va="top", fontsize=11, fontweight="bold", color=INK)
axs[0, 0].set_ylim(0, 290)
fig.legend(handles=[h0, h1], loc="upper right", ncol=2, frameon=False, labelcolor=MUTED)
for ax in axs[1]:
    ax.set_xlabel("Time (s)")
for ax in axs[:, 0]:
    ax.set_ylabel("Inter-story drift (nm)")
fig.suptitle("Damping sensitivity: inter-story drift traces", x=0.01, ha="left", fontsize=14, fontweight="bold")
fig.savefig(FIG / "sensitivity.png", dpi=170)
plt.close(fig)
print("figures ->", FIG)
