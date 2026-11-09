#!/usr/bin/env python3
"""Which nuclear-morphology features differ between malignant and benign breast tumours?

Dataset: Wisconsin Diagnostic Breast Cancer (WDBC), UCI, 569 fine-needle-aspirate images.
Run: python analysis.py        Output goes to results/.
Needs: numpy, pandas, scipy, matplotlib.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ========================== PARAMETERS ==========================
DATA        = "data/wdbc.csv"
OUTDIR      = "results"
ALPHA       = 0.05      # FDR level (Benjamini-Hochberg)
EFFECT_MIN  = 0.80      # AUC threshold for a "large" effect (0.5 = no difference)
N_BOOT      = 2000      # bootstrap resamples for the AUC confidence interval
N_PERM      = 200       # label permutations for the null check
SEED        = 1
# ================================================================

rng = np.random.default_rng(SEED)
df = pd.read_csv(DATA)
y = (df["diagnosis"] == "malignant").to_numpy()
X = df.drop(columns="diagnosis")
print(f"{len(df)} tumours: {y.sum()} malignant, {(~y).sum()} benign; {X.shape[1]} features")

def auc_and_p(x, y):
    """Mann-Whitney U test; AUC = P(malignant value > benign value), ties count half."""
    u, p = stats.mannwhitneyu(x[y], x[~y], alternative="two-sided")
    return u / (y.sum() * (~y).sum()), p

def bh(p):
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)

# 1) one test per feature
rows = []
for c in X.columns:
    x = X[c].to_numpy()
    auc, p = auc_and_p(x, y)
    boots = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if y[i].all() or (~y[i]).all(): continue
        boots.append(auc_and_p(x[i], y[i])[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rows.append((c, np.median(x[y]), np.median(x[~y]), auc, lo, hi, p))
res = pd.DataFrame(rows, columns=["feature", "median_malignant", "median_benign", "AUC", "AUC_lo95", "AUC_hi95", "p"])

# 2) multiple-testing control
m = len(res)
res["q_BH"] = bh(res["p"])
res["p_bonferroni"] = np.minimum(res["p"] * m, 1)
res["direction"] = np.where(res["AUC"] > 0.5, "higher in malignant", "lower in malignant")
res["effect_size_AUC"] = np.maximum(res["AUC"], 1 - res["AUC"])
res["significant_FDR"] = res["q_BH"] < ALPHA
res["large_effect"] = res["significant_FDR"] & (res["effect_size_AUC"] >= EFFECT_MIN)
res = res.sort_values("effect_size_AUC", ascending=False)
res.to_csv(f"{OUTDIR}/feature_tests.tsv", sep="\t", index=False, float_format="%.4g")

# 3) null check: shuffle labels, count false positives at the same FDR threshold
null_hits = []
for _ in range(N_PERM):
    yp = rng.permutation(y)
    pv = [auc_and_p(X[c].to_numpy(), yp)[1] for c in X.columns]
    null_hits.append(int((bh(pv) < ALPHA).sum()))
null_hits = np.array(null_hits)

# 4) how independent are the 30 tests?
corr = X.corr(method="spearman").abs().to_numpy()
iu = np.triu_indices(m, 1)
n_high = int((corr[iu] > 0.9).sum())

summary = [
    f"tumours: {len(df)} ({y.sum()} malignant, {(~y).sum()} benign); features tested: {m}",
    f"significant at FDR {ALPHA} (BH): {int(res.significant_FDR.sum())} / {m}",
    f"significant after Bonferroni at {ALPHA}: {int((res.p_bonferroni < ALPHA).sum())} / {m}",
    f"large effect (FDR-significant and AUC >= {EFFECT_MIN} or <= {1-EFFECT_MIN:.2f}): {int(res.large_effect.sum())}",
    f"label-permutation null ({N_PERM} shuffles): mean FDR-significant features = {null_hits.mean():.2f}, "
    f"max = {null_hits.max()}, shuffles with any hit = {(null_hits > 0).mean():.1%}",
    f"feature pairs with |Spearman rho| > 0.9: {n_high} of {len(iu[0])}",
]
open(f"{OUTDIR}/summary.txt", "w").write("\n".join(summary) + "\n")
print("\n".join(summary))

# 5) figure: effect size with CI, coloured by FDR significance
r = res.sort_values("AUC")
fig, ax = plt.subplots(figsize=(7, 8))
colour = np.where(r.significant_FDR, "#1f77b4", "#9e9e9e")
ax.hlines(range(m), r.AUC_lo95, r.AUC_hi95, color=colour, lw=1.5)
ax.scatter(r.AUC, range(m), c=colour, s=18, zorder=3)
ax.axvline(0.5, color="k", lw=0.8, ls="--")
ax.set_yticks(range(m)); ax.set_yticklabels(r.feature, fontsize=7)
ax.set_xlabel("AUC (0.5 = no difference; > 0.5 = higher in malignant)\nbars: bootstrap 95% CI; blue = significant at FDR 5%")
ax.set_title("Nuclear features: malignant vs benign (WDBC, n = 569)")
fig.tight_layout(); fig.savefig(f"{OUTDIR}/effect_sizes.png", dpi=150)
