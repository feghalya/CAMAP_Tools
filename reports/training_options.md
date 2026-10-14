# CAMAP training options: results summary

## Sources

- **Accuracy**: `output/pytorch-*/<params>-<seed>/accuracy-on-best-validation.json`
  (636 model runs, 12 seeds per setting). Test accuracy on each dataset's own test split.
- **AUC**: `output/revisions/aucs/{homology,binning,shuffling,bias}[.valid].tsv`, computed by
  `utils/aucs.py`. Reported on the BLCL test split and as the mean over the four validation
  cell lines (B721, PBMCs, CT26, EL4).
- **Shuffle types**: `output/stats/auc_shufflings.tsv` (no producing script found).
- **Top of the ranking**: `output/results/final-2021y09m04d/{logregcv,mlp,rf,knn}.ppv.tsv`,
  written by `utils/publish.py` from the final regression
  (`output/metrics/gGRCh37.75gGRCm38.78_pad162_maxRankNP1_maxCont10_maxRedSGD_ANNMethSGDe4000/`),
  which uses the 2020 paper models (`-p e4000`).

Accuracy is only comparable between models evaluated on the same test split. When an option
changes the training dataset itself (context limit, decoy matching), use the validation-set
AUC. Seed-to-seed standard deviation is about 0.004-0.007, so differences below about 0.005
are ties.

## Ranking by option

### 1. Context size (SGD, 2020 paper dataset, test accuracy)

| Rank | Context | Accuracy |
|---|---|---|
| 1 | 243 nt | 0.656 |
| 2 | 162 nt | 0.642 |
| 3 | 81 nt | 0.636 |
| 4 | 27 nt | 0.607 |
| 5 | 9 nt | 0.589 |

The order is the same with shuffled codons and with masking. Gains level off after 162 nt.

### 2. Masking one side (same dataset, test accuracy)

| Context | No mask | Mask after peptide | Mask before peptide |
|---|---|---|---|
| 243 nt | 0.656 | 0.629 | 0.624 |
| 162 nt | 0.642 | 0.618 | 0.613 |

Both sides carry signal, the context after the peptide slightly more at long lengths.

### 3. Real vs. shuffled codons (control)

Real codons always win: about +0.02-0.03 accuracy, and a validation-mean AUC of 0.643 vs. 0.615.

| Shuffle type (BLCL test AUC) | AUC |
|---|---|
| Original | 0.689 |
| Shuffled within region | 0.687 |
| Shuffled 3rd position in region | 0.676 |
| Shuffled genome-wide | 0.647 |

### 4. Optimizer (same datasets, test accuracy)

| Dataset | Adam (500 epochs) | SGD (4000 epochs) |
|---|---|---|
| default (`pep9t5`) | 0.648 | 0.638 |
| no context limit (`pep9t5mINF`) | 0.655 | 0.643 |
| 3 contexts (`pep9t5m3`) | 0.645 | 0.636 |
| TPM-matched (`pep9t5y`) | 0.621 | 0.606 |

Adam beats SGD by about 0.01 on every dataset.

### 5. Context limit (`maxContexts`)

| Limit | BLCL test AUC | Validation mean AUC |
|---|---|---|
| Unlimited (`CINF`) | 0.697 | 0.640 |
| 10 | 0.691 | 0.643 |
| 3 | 0.690 | 0.640 |

Effectively a tie.

### 6. Decoy matching

| Matching | BLCL test AUC | Validation mean AUC |
|---|---|---|
| Binding strength (`sameBS`) | 0.696 | 0.645 |
| None | 0.691 | 0.643 |
| GC content (`sameGC`) | 0.674 | 0.634 |
| Expression (`sameTPM`) | 0.652 | 0.632 |

Matching on binding strength costs nothing; matching on GC or expression lowers the AUC,
since part of what CAMAP captures overlaps with them.

### 7. Bias

| Model | BLCL test AUC | Validation mean AUC | Accuracy |
|---|---|---|---|
| No bias | 0.690 | 0.647 | 0.639 |
| Bias | 0.690 | 0.643 | 0.638 |

A tie. The shuffled no-bias model is the weakest overall (validation mean 0.606).

### 8. Encoding

| Encoding | BLCL test AUC | Validation mean AUC | Accuracy |
|---|---|---|---|
| Codons | 0.691 | 0.643 | 0.638 |
| Amino acids | 0.678 | 0.645 | 0.625 |

Codons win within BLCL; tied on the other cell lines.

### Overall

- **Clear effects:** real codons over shuffled, longer context (243 >= 162 > 81 >> 27, 9),
  no masking, Adam over SGD.
- **Ties:** context limit (3, 10, unlimited), bias vs. no bias, codon vs. amino-acid encoding
  outside BLCL.
- **Decoy matching:** binding strength is neutral; GC or expression costs about 0.01 AUC.

## CAMAP at the top of the ranking

In the final regression, adding CAMAP to binding strength (BS) and expression (TPM) barely
changes the AUC (at most +0.003, and +0.001 over the shuffled-codon control). Its effect is at
the top of the ranking: the false-positive rate among the candidates needed to recover the top
1%, 5% or 10% of true hits (lower is better).

| Dataset (hits at 1%) | Classifier | BS+TPM | + shuffled CAMAP | + CAMAP | CAMAP vs. shuffled at 1% / 5% / 10% |
|---|---|---|---|---|---|
| BLCL (88) | logistic regression | 39.1% | 24.4% | 8.9% | -15.5 / -4.9 / -2.8 points |
| BLCL | MLP | 8.8% | 6.9% | 4.3% | -2.5 / -1.0 / -1.0 |
| BLCL | random forest | 2.3% | 2.2% | 1.1% | -1.1 / -1.0 / -0.5 |
| PBMCs (16) | logistic regression | 88.8% | 85.9% | 83.4% | -2.5 / -2.7 / -1.4 |
| PBMCs | random forest | 52.7% | 48.3% | 36.9% | -11.4 / -5.5 / -3.2 |
| PBMCs | KNN | 59.9% | 51.5% | 41.8% | -9.8 / -4.5 / -2.2 |
| B721 (46) | logistic regression | 52.6% | 45.8% | 41.7% | -4.1 / -2.2 / -2.9 |

- **Human cell lines:** CAMAP lowers false positives at the top with almost every classifier,
  beyond the shuffled control, so the gain is codon-specific. By the top 25-50% the difference
  is about 0 everywhere, which is why the AUC gain is negligible.
- **Part of the gain is not codon-specific:** shuffled CAMAP already removes some false
  positives (BLCL, logistic regression: 39% to 24%, then 9% with real codons).
- **Mouse (EL4, CT26):** inconclusive. The top 1% is only 4-7 hits, and CAMAP helps or hurts
  depending on the classifier.
- **Noise:** the top 1% holds few hits (PBMCs 16, B721 46), so single values are noisy.

### Consequences for choosing training options

- **The rankings above measure the wrong thing for CAMAP's main use.** The metric that matters
  is the false-positive rate (or PPV) at the top 1-5%, on top of BS + TPM, compared with the
  shuffled control, on held-out cell lines.
- **Decoy matching may rank differently on that metric.** `sameTPM` and `sameGC` models score
  lower AUC, but in the regression TPM is already a feature, so only signal independent of
  expression can help at the top. A TPM-matched model could therefore do better there.
- **The experiment that answers it:** rerun step 8 with `-p e4000pep9t5y` (TPM-matched) and
  `-p e4000pep9t5z` (BS-matched) and compare PPV01 and PPV05 with the `-p e4000` results. The
  proteome is already annotated with both model sets (step 7).

## GRCh38 exploration datasets (`output/trainDS.GRCh38.unused/`, now `output/depr-trainDS/`, deleted after this report)

Exploration grid from June-July 2020, built on GRCh38.98 / GRCm38.98 before the final settings
(GRCh37.75 / GRCm38.78, `maxBS500`, `maxContexts10`, `ratio5`, `peplen9`, 12 seeds) were chosen.
354 files, 86 datasets; each dataset has Metadata, CodonEmbeddings and CodonShuffleEmbeddings
files. All use `padding162`. No model, training log or result references them.

### Parameters tested

| Parameter | Values |
|---|---|
| Cell line | BLCL, B721, CT26, PBMC, EL4 |
| Genome | GRCm38.98, GRCh38.98 |
| NetMHC rank threshold (`maxRank`) | 1, 2, 100, 250, 500, 999, 1000, 1250, 2500 |
| Contexts per peptide (`maxContexts`) | 1, 2, 3, 10, 11, 20, 100, 100000 |
| Decoys per hit (`ratio`) | 1, 5, 10, 15, 20 |
| Peptide length (`peplen`) | 8, 9, 10, 11 |
| Filters | `maxTPM`, `minTPM`, `noMaxContextsNS`, `nsMaxTPM`, `nsMinTPM`, `nsTrainMaxTPM`, `nsTrainMinTPM`, `nsTrainValMaxBS`, `nsTrainValMaxTPM`, `nsTrainValMinTPM`, `nsrank`, `nsratio`, `nstpm`, `sameTPM`, `tpmsplit` |

### Combinations (one row per dataset)

| Dataset | Seeds |
|---|---|
| `B721_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_maxRank1_padding162_ratio5_peplen9` | 1-3 |
| `BLCL_GRCh38.98_maxRank1_padding162_ratio5` | 1-12 |
| `BLCL_GRCh38.98_maxRank2_padding162_ratio5` | 1-12 |
| `BLCL_GRCh38.98_padding162_maxRank1000_maxContexts1_ratio10` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1000_maxContexts20_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank100_maxContexts100000_ratio5_peplen9` | 1-3 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts100000_ratio5_peplen9` | 1-3 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts100_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts10_ratio5_noMaxContextsNS_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts1_ratio5_noMaxContextsNS_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts1_ratio5_nsTrainValMaxBS500_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts1_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts1_ratio5` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts20_ratio5_nsTrainMaxTPM0.1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts20_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts2_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio15_sameTPM_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio1_sameTPM_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_maxTPM10_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_maxTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_noMaxContextsNS_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsMaxTPM0.01_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsMaxTPM0.1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsMaxTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsMinTPM0.01_nsMaxTPM0.1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainMaxTPM0.1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainMaxTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainMinTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainValMaxBS50_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainValMaxTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainValMinTPM0.01_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainValMinTPM0.1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainValMinTPM10_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsTrainValMinTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nsratio500_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nstpm10_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nstpm1_nsrank500_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_nstpm1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_sameTPM_peplen9` | 1-3 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_tpmsplit10_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1_maxContexts100000_ratio5_peplen9` | 1-3 |
| `BLCL_GRCh38.98_padding162_maxRank1_maxContexts100_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1_maxContexts1_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank1_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2500_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2500_maxContexts1_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2500_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank250_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank250_maxContexts1_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank250_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2_maxContexts100_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2_maxContexts1_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank2_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts100_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio20_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio20_sameTPM` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio20` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_minTPM0.1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_minTPM1_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_peplen10` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_peplen11` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_peplen8` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_sameTPM_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5_sameTPM` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts10_ratio5` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts11_ratio5_peplen10` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts11_ratio5_peplen11` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts1_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank500_maxContexts3_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank999_maxContexts10_ratio5_peplen9` | 1 |
| `BLCL_GRCh38.98_padding162_maxRank999_maxContexts20_ratio5_peplen9` | 1 |
| `CT26_GRCm38.98_padding162_maxRank1250_maxContexts10_ratio5` | 1 |
| `CT26_GRCm38.98_padding162_maxRank1250_maxContexts1_ratio5_peplen9` | 1 |
| `CT26_GRCm38.98_padding162_maxRank1250_maxContexts1_ratio5` | 1 |
| `EL4_GRCm38.98_padding162_maxRank1250_maxContexts10_ratio5` | 1 |
| `EL4_GRCm38.98_padding162_maxRank1250_maxContexts1_ratio5_peplen9` | 1 |
| `EL4_GRCm38.98_padding162_maxRank1250_maxContexts1_ratio5` | 1 |
| `PBMC_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5_peplen9` | 1 |
| `PBMC_GRCh38.98_padding162_maxRank1250_maxContexts3_ratio5` | 1 |
