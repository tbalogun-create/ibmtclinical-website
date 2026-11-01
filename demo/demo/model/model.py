"""
IBMT Clinical — Pharmacogenomic Dose Simulator: calibration and validation.

Reproduces every number shown on the published page. Run: python3 model.py

Two gene-drug pairs, fitted independently on published concentration data, then
checked against dosing guidelines that were NOT used in the fit.

  Efavirenz   CYP2B6 516G>T x 983T>C (*18)   -> CPIC 2019
  Tacrolimus  CYP3A5 *1/*3                   -> CPIC 2015

Sources
  [1] Sci Rep 2017 (PMC5701031), n=456. Efavirenz C12 by 516 genotype.
  [2] PMC8576856, n=262. CYP2B6 516 genotype prevalence, Jos, Nigeria.
  [3] Bienczak 2016, Br J Clin Pharmacol, n=169. Composite 516x983 clearance ratios.
  [4] PMC2000644, n=170. 983T>C frequency in West Africa; 0% in 93 Caucasian-Americans.
  [5] PMC6739160. CPIC guideline, CYP2B6 and efavirenz.
  [6] s00228-015-1990-0, n=164. Tacrolimus dose-normalised trough by CYP3A5 genotype.
  [7] PMC3655848. CYP3A5*3 frequency across African populations.
  [8] PMC4481158. CPIC guideline, CYP3A5 and tacrolimus.
"""
import numpy as np
from scipy import stats, optimize

HOURS = 24.0
rule = lambda t: print("\n" + t + "\n" + "=" * len(t))

# ============================================================================
# 1. EFAVIRENZ — fit CYP2B6 516 clearances on four published quantities [1]
# ============================================================================
D_GT, D_TT = 0.976, 4.871        # mean C12 above GG at 600 mg, mg/L
P2_GG, P2_GT = 114/306, 16/142   # fraction below 2 mg/L
MEDIAN_ALL = 2.41                # whole-cohort median C12
MIX_FIT = {"GG": 306/456, "GT": 142/456, "TT": 8/456}

ln = np.log
c12_mu    = lambda dose, cl: ln(dose / HOURS) - ln(cl)
mean_c12  = lambda dose, cl, s: np.exp(c12_mu(dose, cl) + s**2 / 2)
frac_below = lambda x, dose, cl, s: stats.norm.cdf((ln(x) - c12_mu(dose, cl)) / s)

def residuals(p):
    cl_gg, cl_gt, cl_tt, s = p
    m = lambda cl: mean_c12(600, cl, s)
    return [ (m(cl_gt) - m(cl_gg)) - D_GT,
             (m(cl_tt) - m(cl_gg)) - D_TT,
             (frac_below(2, 600, cl_gg, s) - P2_GG) * 5,
             (frac_below(2, 600, cl_gt, s) - P2_GT) * 5 ]

CL_GG, CL_GT, CL_TT, SIGMA = optimize.least_squares(
    residuals, x0=[11.0, 7.5, 3.7, 0.45],
    bounds=([4, 3, 1, 0.15], [25, 20, 12, 1.2])).x

rule("EFAVIRENZ — fitted parameters (4 parameters, 4 observations)")
print(f"  CL/F  516GG = {CL_GG:6.3f} L/h")
print(f"  CL/F  516GT = {CL_GT:6.3f} L/h")
print(f"  CL/F  516TT = {CL_TT:6.3f} L/h")
print(f"  lognormal sigma = {SIGMA:.4f}  (CV {np.sqrt(np.exp(SIGMA**2)-1)*100:.0f}%)")

rule("EFAVIRENZ — calibration targets at 600 mg")
m_gg = mean_c12(600, CL_GG, SIGMA)
for name, cl, dtar, ptar in [("GG", CL_GG, 0.0, P2_GG),
                             ("GT", CL_GT, D_GT, P2_GT),
                             ("TT", CL_TT, D_TT, 0.0)]:
    m, p2 = mean_c12(600, cl, SIGMA), frac_below(2, 600, cl, SIGMA)
    p4 = 1 - frac_below(4, 600, cl, SIGMA)
    print(f"  {name}: mean {m:5.2f} (diff {m-m_gg:+5.2f} vs target {dtar:+5.2f}) | "
          f"<2 mg/L {p2*100:5.1f}% vs {ptar*100:4.1f}% | >=4 mg/L {p4*100:5.1f}%")

# ----------------------------------------------------------------------------
# 2. EFAVIRENZ — add 983T>C (*18) via Bienczak composite ratios [3]
#    Their clearances are paediatric, scaled to 15.4 kg, so only the RATIOS
#    transfer. 516-only genotypes keep the adult values fitted above.
# ----------------------------------------------------------------------------
BIEN = {"GG|TT": 6.94, "GT|TT": 4.90, "GG|TC": 3.93,
        "TT|TT": 1.92, "GT|TC": 1.36, "GG|CC": 0.74}
ratio = {k: v / BIEN["GG|TT"] for k, v in BIEN.items()}

EFV_CL = {
    "GG|TT": CL_GG,                       # fitted on adult data [1]
    "GT|TT": CL_GT,                       # fitted
    "TT|TT": CL_TT,                       # fitted
    "GG|TC": CL_GG * ratio["GG|TC"],      # Bienczak ratio x fitted GG [3]
    "GT|TC": CL_GG * ratio["GT|TC"],
    "slow":  CL_GG * ratio["GG|CC"],      # GG|CC and rarer combinations
}

rule("EFAVIRENZ — independent cross-check of the 516 effect")
print("  Fitted on Taiwanese adults [1] vs ratios from Zambian/Ugandan children [3]:")
for g, fitted in [("GT|TT", CL_GT), ("TT|TT", CL_TT)]:
    derived = CL_GG * ratio[g]
    print(f"    {g}: fitted {fitted:5.3f}  Bienczak-derived {derived:5.3f}  "
          f"({abs(derived-fitted)/fitted*100:4.1f}% apart)")
print("  Neither cohort informed the other. The GT agreement is the stronger signal.")

def efv_mix(f516, q983):
    """Composite genotype prevalence, assuming the two loci are independent."""
    tt, tc, cc = (1-q983)**2, 2*q983*(1-q983), q983**2
    m = {"GG|TT": f516["GG"]*tt, "GT|TT": f516["GT"]*tt, "TT|TT": f516["TT"]*tt,
         "GG|TC": f516["GG"]*tc, "GT|TC": f516["GT"]*tc,
         "slow":  f516["TT"]*tc + sum(f516.values())*cc}
    tot = sum(m.values())
    return {k: v/tot for k, v in m.items()}

hw = lambda q: {"GG": (1-q)**2, "GT": 2*q*(1-q), "TT": q*q}

EFV_POPS = {                                     # 516 source, 983 source
    "Nigeria":           (dict(GG=.359, GT=.466, TT=.176), 0.047),   # [2], [4]
    "Uganda & Zambia":   (hw(0.39),                        0.070),   # [3]
    "African American":  (hw(0.33),                        0.075),   # [4]
    "European ancestry": (dict(GG=.563, GT=.375, TT=.063), 0.000),   # gnomAD NFE, [4]
}

def efv_stats(dose, mix):
    cdf = lambda x: sum(w * frac_below(x, dose, EFV_CL[g], SIGMA) for g, w in mix.items())
    lo, hi = cdf(1), cdf(4)
    return dict(sub=lo*100, win=(hi-lo)*100, tox=(1-hi)*100)

rule("EFAVIRENZ — cohorts at 600 mg, with and without 983T>C")
for name, (f516, q983) in EFV_POPS.items():
    off = efv_stats(600, efv_mix(f516, 0.0))
    on  = efv_stats(600, efv_mix(f516, q983))
    note = "unchanged, variant absent" if q983 == 0 else f"{off['tox']:.0f}% -> {on['tox']:.0f}%"
    print(f"  {name:18s} above 4 mg/L: {note}")

rule("EFAVIRENZ — dose by genotype vs CPIC 2019 [5] (not used in fitting)")
AVAILABLE = [200, 400, 600, 800]                 # buildable from 200/600 mg tablets
CPIC_EFV = {"GG|TT": "600 mg", "GT|TT": "400 mg", "GG|TC": "400 mg",
            "TT|TT": "400 or 200 mg", "GT|TC": "400 or 200 mg", "slow": "400 or 200 mg"}
for g in EFV_CL:
    best = max(AVAILABLE, key=lambda d: efv_stats(d, {g: 1.0})["win"])
    ok = "match" if str(best) in CPIC_EFV[g] else "DIFFERS"
    print(f"  {g:7s} -> {best:3d} mg   CPIC: {CPIC_EFV[g]:14s} {ok}")

# ============================================================================
# 3. TACROLIMUS — CYP3A5, calibrated on dose-normalised trough [6]
#    Opposite direction: expressers clear faster and need MORE drug.
# ============================================================================
TAC_CD = {"EXP": 54.61, "NON": 91.80}            # ng/mL per mg/kg/day, day 3 [6]
TAC_IQR = {"EXP": (31.98, 78.87), "NON": (57.60, 130.20)}
TAC_SIG = {g: ln(q3/q1) / (2*0.674490) for g, (q1, q3) in TAC_IQR.items()}
LO, HI, SAFE_LO, SAFE_HI = 10, 15, 4, 20         # protocol NCT04825262 bounds

rule("TACROLIMUS — sigma derived from published interquartile ranges")
for g in TAC_CD:
    s = TAC_SIG[g]
    q3 = np.exp(ln(TAC_CD[g]) + 0.674490*s)
    print(f"  {g:4s} median {TAC_CD[g]:6.2f}  sigma {s:.4f} "
          f"(CV {np.sqrt(np.exp(s*s)-1)*100:.0f}%)  Q3 {q3:6.1f} vs published {TAC_IQR[g][1]}")

def tac_stats(dose, mix):
    cdf = lambda x: sum(w * stats.norm.cdf((ln(x) - ln(dose*TAC_CD[g])) / TAC_SIG[g])
                        for g, w in mix.items())
    lo, hi = cdf(LO), cdf(HI)
    return dict(below=lo*100, win=(hi-lo)*100, above=(1-hi)*100,
                outside=(cdf(SAFE_LO) + 1 - cdf(SAFE_HI))*100)

expresser = lambda q3: {"EXP": 1 - q3*q3, "NON": q3*q3}

TAC_POPS = {                          # CYP3A5*3 frequency, ethnic sample [7]
    "Uganda":            0.04,  "Ghana":   0.11,  "Nigeria": 0.13,
    "Senegal":           0.25,  "Ethiopia": 0.67, "Morocco": 0.80,
    "European ancestry": 0.92,
}
TAC_AA = {"EXP": 0.902, "NON": 0.098}            # observed genotypes, n=357

rule("TACROLIMUS — why 'Africa' is not one population [7]")
for name, q3 in TAC_POPS.items():
    print(f"  {name:18s} CYP3A5*3 {q3:.2f}  ->  {(1-q3*q3)*100:5.1f}% expressers")
print(f"  {'African American':18s} CYP3A5*3 0.29  ->  {TAC_AA['EXP']*100:5.1f}% expressers")
print("  Twentyfold spread within Africa. Ethiopia sits closer to Europe than to Nigeria.")

rule("TACROLIMUS — Nigeria at the protocol's non-expresser dose (0.15 mg/kg/day)")
ng, eu = expresser(0.13), expresser(0.92)
for name, mix in [("Nigeria", ng), ("European ancestry", eu)]:
    st = tac_stats(0.15, mix)
    print(f"  {name:18s} below target {st['below']:5.1f}% | in target {st['win']:5.1f}% | "
          f"outside {SAFE_LO}-{SAFE_HI} {st['outside']:5.1f}%")

rule("TACROLIMUS — dose vs CPIC 2015 [8] (not used in fitting)")
grid = [round(d, 2) for d in np.arange(0.05, 0.401, 0.01)]
free = max(grid, key=lambda d: tac_stats(d, ng)["win"])
safe = max((d for d in grid if tac_stats(d, ng)["outside"] <= 25),
           key=lambda d: tac_stats(d, ng)["win"])
non  = max(grid, key=lambda d: tac_stats(d, {"NON": 1.0})["win"])
exp_ = max(grid, key=lambda d: tac_stats(d, {"EXP": 1.0})["win"])
print(f"  Expresser optimum {exp_:.2f} / non-expresser optimum {non:.2f} "
      f"= {exp_/non:.2f}x   CPIC says 1.5-2x  -> "
      f"{'match' if 1.5 <= exp_/non <= 2 else 'DIFFERS'}")
print(f"  Nigeria, unconstrained: {free:.2f} mg/kg/day, "
      f"{tac_stats(free, ng)['outside']:.1f}% outside bounds")
print(f"  Nigeria, respecting the protocol's 25% stopping rule: {safe:.2f} mg/kg/day, "
      f"{tac_stats(safe, ng)['outside']:.1f}% outside")
print(f"  The stopping rule caps the dose at {safe/0.15:.2f}x the non-expresser dose, "
      f"BELOW the 1.5-2x CPIC range.")
print("  Finding: in a ~98% expresser cohort this protocol's own safety bound")
print("  prevents the starting dose that clinical guidance calls for.")

rule("Constants used by the published page")
print(f"  EFV_CL = {{" + ", ".join(f'"{k}":{v:.4f}' for k, v in EFV_CL.items()) + "}")
print(f"  EFV_SIGMA = {SIGMA:.4f}")
print(f"  TAC_CD = {TAC_CD}")
print(f"  TAC_SIG = {{" + ", ".join(f'"{k}":{v:.4f}' for k, v in TAC_SIG.items()) + "}")
