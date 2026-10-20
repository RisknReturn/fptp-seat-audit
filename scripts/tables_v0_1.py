"""Markdown tables for the report, generated from out/*.json (so the doc numbers come from the runs)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import json, numpy as np, sys
sys.path.insert(0, _ROOT)
from ukseat import config_2024 as C
OUT = f"{_ROOT}/out"
bt = json.load(open(f"{OUT}/backtest2_v0_1.json"))
cv = json.load(open(f"{OUT}/converge_v0_1.json"))
mc = json.load(open(f"{OUT}/mc_compare_std_v0_1.json"))
ref = json.load(open(f"{OUT}/ref_exact2_v0_1.json"))
ref_old = json.load(open(f"{OUT}/ref_exact_v0_1.json"))
l2 = json.load(open(f"{OUT}/nat_ladder2_v0_1.json"))
ACT = C.ACTUAL_SEATS_GB_2024
P = ["Lab", "Con", "LD", "Ref", "SNP", "Grn", "PC", "Oth"]
L = []
def p(s=""): L.append(s)

p(f"### T1. 2024 백테스트: 정당별 의석 (λ=0, 기준 격자: 전국 {ref['nodes']}노드 × 지역 {ref['reg_nodes'][0]}/{ref['reg_nodes'][5]}노드)")
p("| 정당 | 실제 | 평균 | 표준편차 | 5% | 50% | 95% | P(S ≤ 실제) |")
p("|---|---:|---:|---:|---:|---:|---:|---:|")
dist = np.array(ref["dist"])
for k in P:
    s = ref["summ"][k]; i = ["Con","Lab","LD","Ref","Grn","SNP","PC","Oth"].index(k)
    pit = dist[i, : ACT[k] + 1].sum()
    p(f"| {k} | {ACT[k]} | {s['mean']:.1f} | {s['sd']:.1f} | {s['q05']} | {s['q50']} | {s['q95']} | {pit:.3f} |")
p()
p("### T2. 주요 사건 확률 (λ=0)")
p("| 사건 | 정확 계산 (기준) | MC 3×10^6 | (정확−MC)/SE | 이전 기준 격자(1차, 축 배분 오류) | (이전−MC)/SE |")
p("|---|---:|---:|---:|---:|---:|")
m6 = [r for r in mc["runs"] if r["N"] == 1_000_000]
for k, v in ref["ev"].items():
    mv = np.mean([r["ev"][k] for r in m6])
    se = np.sqrt(max(mv * (1 - mv), 1e-12) / 3e6)
    vo = ref_old["ev"][k]
    p(f"| {k} | {v:.5f} | {mv:.5f} | {(v-mv)/se:+.1f} | {vo:.5f} | {(vo-mv)/se:+.1f} |")
p()
p("### T3. 구적 격자 사다리: (a) 전국 격자(지역 R0 고정, 기준 N6×R0), (b) 지역 격자(전국 N2 고정, 기준 N2×R5)")
p("| 사다리 | 격자 | 전국 노드 | 지역 노드(잉글랜드 한 지역/스코틀랜드) | 시간(초, 2코어) | 최대 |오차| (18개 사건) | 그 사건 |")
p("|---|---|---:|---:|---:|---:|---|")
base = {"nat": next(r for r in cv if r["ladder"] == "nat" and r["nat"] == "N6"),
        "reg": next(r for r in cv if r["ladder"] == "reg" and r["reg"] == "R5")}
for r in cv:
    b0 = base[r["ladder"]]
    d = {k: abs(r["ev"][k] - b0["ev"][k]) for k in b0["ev"]}
    km = max(d, key=d.get)
    p(f"| {'전국' if r['ladder']=='nat' else '지역'} | {r['nat']}×{r['reg']} | {r['nodes']} | {r['reg_nodes'][0]}/{r['reg_nodes'][5]} | {r['sec']:.1f} | {d[km]:.1e} | {km} |")
p(f"1차 사다리의 한계: 전국 4·5번 축과 지역 4~6번 축을 세분하지 않아, 이 표로 추정한 오차(최대 6e-4)는 실제 오차(SNP≤9에서 0.014, Oth≥5에서 0.016, T2 '이전 기준' 열)를 크게 과소평가했습니다.")
p(f"현재 기준 격자: 전국(정당별 정규화 회전) 사다리꼴 h=0.35 × GH 7/5/5/3 = {ref['nodes']}노드, 지역 GH (5,5,3,2,2,2) = 잉글랜드 지역당 {ref['reg_nodes'][0]} / 스코틀랜드·웨일스 {ref['reg_nodes'][5]}노드, {ref['sec']:.0f}초. 절단 오차 상한(TV) {ref['tv']:.1e}.")
p()
p("### T3b. 전국 축을 하나씩 세분 (지역 R1 고정, MC 3×10^6 대비 최대 |차이|)")
p("| 회전 기준 | 격자 | 전국 노드 | 시간(초) | 최대 |정확−MC| | 그 사건 |")
p("|---|---|---:|---:|---:|---|")
mce = l2["mc"]
for r in l2["rows"]:
    d = {k: abs(r["ev"][k] - mce[k]) for k in mce}; km = max(d, key=d.get)
    p(f"| {'정당별 정규화 영향' if r['norm'] else '의석 영향(원값)'} | {r['tag']} | {r['nodes']} | {r['sec']:.0f} | {d[km]:.4f} | {km} |")
p()
p("### T4. 몬테카를로 (같은 모델의 생성 과정을 직접 시뮬레이션)")
p("| N | 반복 | 평균 시간(초) | RMSE P(Lab≥400) | RMSE P(Lab<326) | 상대 RMSE P(Lab≤250) | 상대 RMSE P(Con≥250) |")
p("|---:|---:|---:|---:|---:|---:|---:|")
for N in sorted(set(r["N"] for r in mc["runs"])):
    rr = [r for r in mc["runs"] if r["N"] == N]
    f = lambda ev: np.sqrt(np.mean([(r["ev"][ev] - ref["ev"][ev]) ** 2 for r in rr]))
    p(f"| {N:,} | {len(rr)} | {np.mean([r['sec'] for r in rr]):.2f} | {f('Lab>=400'):.1e} | {f('Lab>=326 (majority)'):.1e} | "
      f"{f('Lab<=250')/ref['ev']['Lab<=250']:.2f} | {f('Con>=250')/ref['ev']['Con>=250']:.2f} |")
p()
p("### T5. 시나리오 차이 λ=0.5 − λ=0.25: MC(공통 난수, N=10^5, 10회)의 흩어짐 vs 정확 계산")
ex = {k: bt["lam0.5"]["events"][k] - bt["lam0.25"]["events"][k] for k in bt["lam0.5"]["events"]}
p("| 사건 | 정확 차이 | MC 차이 평균 | MC 차이 표준편차 |")
p("|---|---:|---:|---:|")
for k in ["Lab>=400", "Lab>=450", "Lab<=411 (actual or fewer)", "Con>=121 (actual or more)", "Con<100", "LD>=72 (actual or more)", "Hung (no party>=326)"]:
    xs = [d[k] for d in mc["crn_diffs_N1e5"]]
    p(f"| {k} | {ex[k]:+.4f} | {np.mean(xs):+.4f} | {np.std(xs, ddof=1):.4f} |")
p()
p("### T6. λ 스캔 (2015년 크기 여론조사 편향의 배수; 격자: 정규화 회전, 사다리꼴 0.35 × GH 5/3/5/3, 지역 (3,3,2,2,2,2))")
p("| λ | 중심 GB 득표율 Con/Lab/LD/Ref/Grn | E[Lab] | E[Con] | E[LD] | E[Ref] | P(Lab≤411) | P(Con≥121) | P(Lab≥326) |")
p("|---:|---|---:|---:|---:|---:|---:|---:|---:|")
for l in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]:
    r = bt[f"lam{l}"]; g = r["gb_shares_centre"]; s = r["summary"]; e = r["events"]
    p(f"| {l} | {g['Con']:.1f}/{g['Lab']:.1f}/{g['LD']:.1f}/{g['Ref']:.1f}/{g['Grn']:.1f} | {s['Lab']['mean']:.0f} | {s['Con']['mean']:.0f} | "
      f"{s['LD']['mean']:.0f} | {s['Ref']['mean']:.0f} | {e['Lab<=411 (actual or fewer)']:.3f} | {e['Con>=121 (actual or more)']:.3f} | {e['Lab>=326 (majority)']:.3f} |")
p()
p("### T7. 오라클(실제 득표율 입력) vs 최종 여론조사 입력")
p("| 정당 | 실제 | 여론조사 입력 E[S] (5~95%) | 실제 득표율 입력·전국 오차 끔 E[S] (5~95%) |")
p("|---|---:|---|---|")
for k in P:
    a = bt["lam0_ref"]["summary"][k]; o = bt["oracle_no_nat"]["summary"][k]
    p(f"| {k} | {ACT[k]} | {a['mean']:.0f} ({a['q05']}~{a['q95']}) | {o['mean']:.0f} ({o['q05']}~{o['q95']}) |")
p()
p("### T8. 선거구 단위 보정도")
for tag, nm in (("calib_lam0", "여론조사 입력"), ("calib_oracle", "실제 득표율 입력")):
    c = bt[tag]
    p(f"- {nm}: Brier {c['brier']:.3f}, 로그 점수 {c['logscore']:.3f}, 최유력 후보 적중 {c['hit_rate']*631:.0f}/631 (모델 기대 {c['expected_hits']:.0f})")
p()
p("| 예측 확률 구간 | 쌍 수 | 평균 예측 | 실제 비율 |  (여론조사 입력, 선거구×정당 모든 쌍)")
p("|---|---:|---:|---:|")
for r in bt["calib_lam0"]["reliability_all"]:
    p(f"| {r['lo']:.2f}~{r['hi']:.2f} | {r['n']} | {r['mean_p']:.3f} | {r['freq']:.3f} |")
p()
p("가장 의외였던 결과 10곳 (여론조사 입력 기준, 실제 당선자에게 준 확률 낮은 순)")
p("| 선거구 | 지역 | 모델 최유력 (확률) | 실제 당선 | 실제 당선자 확률 |")
p("|---|---|---|---|---:|")
for u in bt["calib_lam0"]["upsets"][:10]:
    p(f"| {u['seat']} | {u['region']} | {u['fav']} ({u['p_fav']:.2f}) | {u['winner']} | {u['p_winner']:.4f} |")
p()
p("### T9. 분산 분해 (λ=0): 전국 / 지역 / 선거구")
p("| 정당 | 전국 | 지역 | 선거구 |")
p("|---|---:|---:|---:|")
for k in P:
    v = bt["lam0_ref"]["var_share"][k]
    p(f"| {k} | {v['nat']:.2f} | {v['reg']:.2f} | {v['seat']:.2f} |")
open(f"{OUT}/tables_v0_1.md", "w").write("\n".join(L))
print("\n".join(L))
