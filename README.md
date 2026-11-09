# Nuclear morphology of malignant vs benign breast tumours (WDBC)

## Biological question
Fine-needle-aspirate (FNA) images of breast masses were digitised and the cell nuclei measured.
**Which nuclear features differ between malignant and benign tumours, and how large are the differences?**
This is a descriptive, exploratory question. It is not a diagnostic model.

## Data
Wisconsin Diagnostic Breast Cancer (WDBC) dataset, UCI Machine Learning Repository
(Street, Wolberg & Mangasarian, 1993): 569 tumours (212 malignant, 357 benign), 30 numeric features.
The 30 features are 10 nuclear measurements (radius, texture, perimeter, area, smoothness, compactness,
concavity, concave points, symmetry, fractal dimension), each summarised three ways per image:
mean, standard error ("error") and worst (mean of the three largest values).
The file `data/wdbc.csv` is the copy distributed with scikit-learn, saved with readable column names.

## How to run
```bash
pip install -r requirements.txt
python analysis.py
```
Takes about a minute. The random seed is fixed, so the output is reproducible. All parameters (FDR level, effect-size threshold, resamples, seed) are at the top of `analysis.py`.

## Methods
1. **Per-feature test:** two-sided Mann-Whitney U test, malignant vs benign (no normality assumption).
2. **Effect size:** AUC, the probability that a random malignant tumour has a higher value than a random benign one (0.5 = no difference). 95% confidence interval from 2,000 bootstrap resamples.
3. **Multiple testing:** 30 tests, so p-values are adjusted with Benjamini-Hochberg (FDR 5%). Bonferroni is reported for comparison.
4. **Sanity check of the error control:** the diagnosis labels were shuffled 200 times and the same procedure was repeated, to see how many "significant" features appear when no real difference exists.
5. **"Large effect"** is defined as FDR-significant and AUC >= 0.80 (or <= 0.20). The 0.80 threshold is my own convention, set in advance and not tuned.

## Results
(`results/feature_tests.tsv`, `results/summary.txt`, `results/effect_sizes.png`)

- **27 of 30 features** differ at FDR 5% (26 after Bonferroni). **15** also have a large effect.
- The strongest differences are in **nucleus size and shape irregularity**: worst perimeter (AUC 0.976), worst radius (0.970), worst area (0.970), worst concave points (0.967), mean concave points (0.964). Malignant nuclei are larger and have more concave points and more concavity. For example, median worst perimeter is 138 in malignant and 87 in benign tumours.
- **Texture, smoothness, symmetry and fractal dimension** differ significantly but modestly (AUC about 0.62 to 0.78).
- **No evidence of a difference** for texture error, mean fractal dimension and smoothness error (q > 0.2).
- **Symmetry error** is only marginally significant (q = 0.031, AUC 0.44, confidence interval reaching 0.49). I would not read anything into it.
- **Null check:** with shuffled labels, the average number of "significant" features was 0.44, and 6.5% of shuffles produced any hit, close to the nominal 5%. When a shuffle did produce a false hit it often produced many (maximum 20), because the features are strongly correlated (see below). The FDR procedure is therefore working, but false positives arrive in clusters, not one by one.

## Interpretation
Malignant tumours in this dataset have nuclei that are larger and more irregular in outline: bigger, with more and deeper concavities. This agrees with the long-standing pathology observation that cancer cell nuclei are enlarged and pleomorphic, which is why the dataset was built. Radius, perimeter and area measure essentially the same geometric property, so they count as **one finding (nuclear size), not three**.

## Limitations and confounders
- **Dependent tests.** 26 of 435 feature pairs have |Spearman rho| > 0.9. The effective number of independent tests is much smaller than 30, so "27 significant features" overstates how many independent signals there are. BH remains valid under positive dependence, but the count should not be read as 27 separate discoveries.
- **P-values are tiny because n is large.** With 569 samples almost any difference is "significant". The effect sizes and confidence intervals carry the meaning, not the q-values.
- **Possible non-independence of samples.** The public file has no patient identifiers. If some patients contributed more than one tumour or image, the effective sample size is smaller and the p-values are too optimistic.
- **Single source, old data.** One institution (University of Wisconsin), collected in the early 1990s. Differences in staining, imaging equipment and population could change the numbers elsewhere. There is no external validation.
- **Confounders not recorded.** Patient age, tumour stage and grade, and sample preparation are unavailable. Tumour size is related to stage, so some of the size signal may reflect more advanced disease and not malignancy itself.
- **Label and measurement issues.** The nuclei were outlined semi-automatically, so segmentation choices may affect the values. The benign/malignant label comes from later diagnosis, with its own error rate.
- **Association, not causation, and not a diagnostic tool.** AUC values are per feature on the same data used to select them. They show that the features are informative, not how well a classifier would perform on new patients. This analysis should not be used for any clinical decision.

## Files
| File | Content |
|---|---|
| `analysis.py` | the whole analysis |
| `data/wdbc.csv` | dataset |
| `requirements.txt` | Python packages |
| `results/feature_tests.tsv` | per-feature medians, AUC, CI, p, q (BH), Bonferroni p |
| `results/summary.txt` | headline numbers and null check |
| `results/effect_sizes.png` | effect sizes with confidence intervals |

## Reference
W.N. Street, W.H. Wolberg, O.L. Mangasarian. Nuclear feature extraction for breast tumor diagnosis. IS&T/SPIE 1993 International Symposium on Electronic Imaging, 1905, 861-870.
