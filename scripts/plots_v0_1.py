"""Figures for the v0.1 report (static PNG; one series per panel where colour would otherwise carry party)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, _ROOT)
from ukseat.model import CATS
from ukseat import config_2024 as C

OUT = f"{_ROOT}/out"; FIG = f"{_ROOT}/figs"
_os.makedirs(FIG, exist_ok=True)
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3 (reference palette)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150,
                     "font.family": ["Noto Sans CJK KR", "DejaVu Sans"]})
bt = json.load(open(f"{OUT}/backtest2_v0_1.json"))
ACT = C.ACTUAL_SEATS_GB_2024
PARTIES = ["Lab", "Con", "LD", "Ref", "SNP", "Grn"]
XL = {"Lab": (200, 560), "Con": (0, 320), "LD": (0, 140), "Ref": (0, 160), "SNP": (0, 57), "Grn": (0, 12)}

# ---- Fig 1: seat distributions, lambda=0 vs lambda=1, actual marked
fig, axs = plt.subplots(2, 3, figsize=(9.5, 5.2))
for ax, p in zip(axs.ravel(), PARTIES):
    i = CATS.index(p)
    for key, col, lab in (("lam0_ref", S1, "최종 여론조사 그대로 (λ=0)"), ("lam1.0", S2, "2015 크기 편향 보정 (λ=1)")):
        d = np.array(bt[key]["dist"])[i]
        ax.plot(np.arange(len(d)), d, color=col, lw=2, label=lab)
    ax.axvline(ACT[p], color=INK, lw=1.2, ls="--")
    ax.text(ACT[p], ax.get_ylim()[1] * 0.92, f" 실제 {ACT[p]}", color=INK, fontsize=8)
    ax.set_xlim(*XL[p]); ax.set_title(p, color=INK, fontsize=10, loc="left"); ax.set_yticklabels([])
h, l = axs[0, 0].get_legend_handles_labels()
fig.suptitle("2024 총선 백테스트: 정당별 의석 분포 (GB 631석, 표본 추출 없는 계산)", color=INK, x=0.01, ha="left")
fig.supxlabel("의석 수 (점선 = 실제)", color=INK2); fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=8.5)
fig.savefig(f"{FIG}/fig1_seat_dist_v0_1.png"); plt.close(fig)

# ---- Fig 2: lambda scan
lams = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
fig, axs = plt.subplots(1, 4, figsize=(10, 2.9))
for ax, p in zip(axs, ["Lab", "Con", "LD", "Ref"]):
    mu = [bt[f"lam{l}"]["summary"][p]["mean"] for l in lams]
    q05 = [bt[f"lam{l}"]["summary"][p]["q05"] for l in lams]; q95 = [bt[f"lam{l}"]["summary"][p]["q95"] for l in lams]
    ax.fill_between(lams, q05, q95, color=S1, alpha=0.18, lw=0)
    ax.plot(lams, mu, color=S1, lw=2, marker="o", ms=4)
    ax.axhline(ACT[p], color=INK, lw=1.2, ls="--"); ax.set_title(f"{p} (점선=실제 {ACT[p]})", loc="left", fontsize=9, color=INK)
    ax.set_xlabel("λ (2015 편향의 배수)")
axs[0].set_ylabel("의석 (평균, 5~95%)")
fig.tight_layout(); fig.savefig(f"{FIG}/fig2_lambda_scan_v0_1.png"); plt.close(fig)

# ---- Fig 3: accuracy vs time (MC replicates vs exact grids), errors relative to the refined reference
try:
    mc = json.load(open(f"{OUT}/mc_compare_std_v0_1.json"))
    ref = json.load(open(f"{OUT}/ref_exact2_v0_1.json"))["ev"]
    cv = json.load(open(f"{OUT}/converge_v0_1.json")); l2 = json.load(open(f"{OUT}/nat_ladder2_v0_1.json"))["rows"]
    EVS = [("Lab>=400", "중심 확률 P(Lab≥400)≈0.48"), ("Lab>=326 (majority)", "꼬리 P(Lab<326)≈0.045"), ("Lab<=250", "깊은 꼬리 P(Lab≤250)≈3e-4")]
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.2))
    for ax, (ev, lab) in zip(axs, EVS):
        Ns = sorted(set(r["N"] for r in mc["runs"]))
        t = [np.mean([r["sec"] for r in mc["runs"] if r["N"] == N]) for N in Ns]
        rm = [np.sqrt(np.mean([(r["ev"][ev] - ref[ev]) ** 2 for r in mc["runs"] if r["N"] == N])) for N in Ns]
        ax.plot(t, rm, color=S2, lw=2, marker="o", ms=5, label="몬테카를로 RMSE")
        nat = [r for r in cv if r["ladder"] == "nat"]
        ax.plot([r["sec"] for r in nat], [max(abs(r["ev"][ev] - ref[ev]), 1e-8) for r in nat], color=S1, lw=2, marker="s", ms=5,
                label="정확 계산, 1차 격자 배분")
        nn = sorted([r for r in l2 if r["norm"]], key=lambda r: r["sec"])
        ax.plot([r["sec"] for r in nn], [max(abs(r["ev"][ev] - ref[ev]), 1e-8) for r in nn], color=S3, lw=0, marker="D", ms=5,
                label="정확 계산, 정규화 배분")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(lab, loc="left", fontsize=9, color=INK)
        ax.set_xlabel("계산 시간 (초, 2코어)")
    axs[0].set_ylabel("기준값 대비 오차 (절대)"); axs[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig3_accuracy_time_v0_1.png"); plt.close(fig)
except FileNotFoundError as e:
    print("skip fig3", e)

# ---- Fig 4: seat-level reliability (all seat-party pairs)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))
for ax, tag, title in ((axs[0], "calib_lam0", "최종 여론조사 입력"), (axs[1], "calib_oracle", "실제 전국 득표율 입력")):
    rel = bt[tag]["reliability_all"]
    x = [r["mean_p"] for r in rel]; y = [r["freq"] for r in rel]; n = [r["n"] for r in rel]
    ax.plot([0, 1], [0, 1], color=INK2, lw=1, ls=":")
    ax.plot(x, y, color=S1, lw=2, marker="o", ms=5)
    ax.set_xlabel("예측 당선 확률"); ax.set_title(f"{title} (Brier {bt[tag]['brier']:.3f})", loc="left", fontsize=9, color=INK)
axs[0].set_ylabel("실제 당선 비율")
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_reliability_v0_1.png"); plt.close(fig)

# ---- Fig 5: variance decomposition
fig, ax = plt.subplots(figsize=(6.4, 2.8))
ps = ["Lab", "Con", "LD", "Ref", "SNP", "Grn"]
vs = bt["lam0_ref"]["var_share"]
left = np.zeros(len(ps))
for comp, col, lab in (("nat", S1, "전국 여론조사 오차"), ("reg", S2, "지역 편차"), ("seat", S3, "선거구 고유")):
    v = np.array([vs[p][comp] for p in ps])
    ax.barh(ps, v, left=left, color=col, edgecolor="white", linewidth=2, label=lab); left += v
ax.invert_yaxis(); ax.set_xlim(0, 1); ax.set_xlabel("의석 수 분산 비중"); ax.legend(frameon=False, ncol=3, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(); fig.savefig(f"{FIG}/fig5_var_decomp_v0_1.png"); plt.close(fig)
print("figs done")
