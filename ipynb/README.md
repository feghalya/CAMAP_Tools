# Notebooks

An archive of the exploratory notebooks behind the CAMAP analyses, collected here from three
projects. They are kept for their saved figures and printed results.

The three `graphs` notebooks are three different analyses, oldest first: `graphs.ipynb`
(Nov 2019, ImmPred), then `Graphs_old.ipynb` and `Graphs_new.ipynb` (Dec 2019,
Codon_publication).

| Notebook | Date | From | What it does | Expects |
|---|---|---|---|---|
| `graphs.ipynb` | 2019-11 | ImmPred | Earliest score-distribution plots, source vs. non-source, with a Mann-Whitney test. 8 figures. | `scores.{source,nonsource}.graham.txt`, `scores.PBMC.*.txt` |
| `Graphs_old.ipynb` | 2019-12 | Codon_publication | AUC analysis per dataset, peptide length and NetMHC version. 3 figures. | `Codon_publication/output/finalTables[.old]/*.tsv.gz` |
| `Graphs_new.ipynb` | 2019-12 | Codon_publication | Reworked AUC analysis, original vs. shuffled sequences by peptide size. 6 figures. | `all_metrics.tsv`, `output/various/metrics.tsv` |
| `camap_stats.ipynb` | 2021-04 | CAMAP_Tools | ROC/AUC and logistic regression, plus the gene-level `max()` aggregation that became `max_reduce`. 3 figures. | `finalTables/*.tsv.gz`, `parsed/peptides.*.all.txt` |
| `camap_regressions.ipynb` | 2021-04 | CAMAP_Tools | Regression comparisons across feature combinations. 9 figures. | `finalTables/*.tsv.gz` |
| `camap.ipynb` | 2021-04 | CAMAP_Tools | **Prototype of the PPV analysis**: how many candidates must be scanned to recover 1%, 10%, ... 99% of true hits, for each feature combination. No figures; the results are in its printed output. | `Human_BLCL.scores.tsv.gz` |

## Why these matter

`camap.ipynb` is the origin of the top-of-ranking analysis. It predates `utils/ppv.py` by four
months and uses the same proportions that became `PPV01`...`PPV99`. It is the earliest record
of the finding that CAMAP's contribution is concentrated in the top few percent of candidates
rather than in the overall AUC. See `reports/training_options.md`.

`camap_stats.ipynb` holds an early version of the per-gene reduction later adopted as
`max_reduce=True` in `Dataset.construct_datasets()`. The mature prototype of that choice,
including the sum-versus-max comparison, is in `utils/prototype.py`.

The rest are superseded: the 2019 AUC work was replaced by `revisions/aucs/`, and the 2021
regression comparisons by `output/metrics/`.
