# TASC user guide: inputs, outputs and interpretation

> **DRAFT** – sections marked `[TO COMPLETE]` still need the authors' input.

1. [Preparing the input](#1-preparing-the-input)
2. [Online app or local run?](#2-online-app-or-local-run)
3. [Understanding the output](#3-understanding-the-output)
4. [Unclassified (NA) and multiple calls: rare subtypes](#4-unclassified-na-and-multiple-calls-rare-subtypes)
5. [The TASC risk score](#5-the-tasc-risk-score)

---

## 1. Preparing the input

TASC takes **raw gene-level counts** (not TPM, FPKM or already-normalised values). It needs one column of gene identifiers in **Ensembl format** (`ENSG…`) and one column per sample. Version suffixes (`ENSG00000227232.5`) are removed automatically. Choose the matching option in *Source of count matrix* before uploading.

### featureCounts (`Source: featureCounts`)
The first six columns are annotation and are skipped. Samples start from column 7.
```
Geneid           Chr   Start  End    Strand  Length  Sample_A  Sample_B
ENSG00000223972  chr1  11869  14409  +       2541    2         0
ENSG00000227232  chr1  14404  29570  -       15167   100       189
```

### Salmon, e.g. nf-core/rnaseq `salmon.merged.gene_counts.tsv` (`Source: Salmon`)
The first two columns (gene ID and gene name) are skipped. Samples start from column 3.
```
gene_id             gene_name  Sample_A  Sample_B
ENSG00000223972.5   DDX11L1    2.000     0.000
ENSG00000227232.5   WASH7P     100.000   189.000
```

### Custom format (`Source: Custom format`)
Gene IDs go in the first column and every other column is a sample.
```
GeneID           Sample_A  Sample_B  Sample_C
ENSG00000139618  523       610       489
ENSG00000157764  342       298       410
```

**Notes**
- The file must be tab-separated (`.tsv` or `.txt`).
- Model genes that are missing from your file are **imputed** from the reference cohort. If more than 10% of model genes are imputed, or if any of the 30 most important genes are imputed, the app shows a warning. Imputation can lower subtype probabilities (see [Section 4](#4-unclassified-na-and-multiple-calls-rare-subtypes)).
- Imputation uses a fixed random seed, so the same file always gives the same result.
- Normalisation (DESeq2 VST) is computed across all the samples in one upload. Upload samples that were processed with the same pipeline together, and avoid mixing very different library types in one file.

## 2. Online app or local run?

The [online app](https://01971d1c-d81d-a598-4238-5afb7b3e381a.share.connect.posit.cloud/) accepts files up to 200 MB and runs on shared resources.

> ⚠️ **If you have more than 100 samples, run TASC on your own computer.** Large uploads in the online app are slow and can time out, and the PCA and metagene plots become hard to read.

To run it locally:
```r
# in a terminal: git clone --depth 1 https://github.com/CBenetti/TASC.git ./TASC
shiny::runApp("path_to_TASC_repo/TASC.R")
```
`TASC.R` is identical to the online app, except that the upload limit is raised to 5 GB.

## 3. Understanding the output

### Predicted Classes table
| Field | Meaning |
|---|---|
| `Sample` | Sample name from the header of your count matrix. |
| `Prediction` | The subtype with the highest probability. One override applies: if the top class is ETP-like and the LMO2 γδ-like probability is above 0.286, the sample is reported as **LMO2 γδ-like**. |

### TASC Risk Score table
| Field | Meaning |
|---|---|
| `Sample` | Sample name. |
| `Prediction` | Same as above. |
| `TASC risk score` | Prognostic score for overall survival. Higher values mean higher risk ([Section 5](#5-the-tasc-risk-score)). |
| `Risk class` | `Low`, `Intermediate` or `High`, based on the tertile cut-points of the reference cohort (−0.0527 and 0.4464). |

### Complete table of probabilities (download)
| Field | Meaning |
|---|---|
| *(13 subtype columns)* | Random Forest vote fraction for each subtype: TLX1, NKX2-1, TLX3, ETP-like, TAL1 DP-like, TAL1 αβ-like, TME-enriched, LMO2 γδ-like, KMT2A, SPI1, MLLT10, HOXA9 TCR, BCL11B. Each row sums to 1. These are relative support scores, not calibrated probabilities. |
| `class` | Subtype with the highest vote fraction, **without** the LMO2 override. |
| `probability` | Vote fraction of `class`. |
| `corrected` | **Thresholded call:** every subtype whose vote fraction is above its subtype-specific threshold (sensitivity and specificity >90%, maximal specificity; Supplementary Table 2 of the paper). `NA` means no subtype passes its threshold. Several subtypes are separated by `;`. |
| `satisfies_prob` | `TRUE` if the top class (`class`) also passes its own threshold. |

Subtype thresholds used for `corrected`:

| Subtype | Threshold | Subtype | Threshold |
|---|---|---|---|
| BCL11B | 0.368 | NKX2-1 | 0.182 |
| ETP-like | 0.377 | SPI1 | 0.386 |
| HOXA9 TCR | 0.188 | TAL1 DP-like | 0.350 |
| KMT2A | 0.300 | TAL1 αβ-like | 0.429 |
| LMO2 γδ-like | 0.089 | TLX1 | 0.235 |
| MLLT10 | 0.327 | TLX3 | 0.336 |
| | | TME-enriched | 0.264 |

**How to read them together:** `Prediction` always names a subtype. `corrected` tells you whether that call is **confident**. A sample whose `Prediction` is ETP-like but whose `corrected` is `NA` should be reported as *unclassified (closest: ETP-like)*.

### Other outputs
- **Imputed gene list:** model genes that were missing from your file and were imputed.
- **PCA plot:** samples in your upload, coloured by `Prediction`.
- **Metagene feature plots:** mean expression of each predicted subtype's signature genes. Signature genes are **not** imputed, so a plot is skipped if too few signature genes are present. *Explore Signature on Demand* lets you check any subtype's signature, which is useful for NA samples.

## 4. Unclassified (NA) and multiple calls: rare subtypes

TASC was trained on the 13 subtypes with at least 10 samples in the reference cohort. Rarer subtypes (**NKX2-5, NUP98, NUP214, STAG2::LMO2** and others) are not model classes. For these samples, `Prediction` is always the *closest* trained subtype, so use the `corrected` field to judge the call.

What TASC reports for the 27 rare-subtype samples of the reference cohort (published thresholds):

| True subtype | n | `corrected` call | Typical closest `Prediction` |
|---|---|---|---|
| STAG2::LMO2 | 8 | NA 8/8 | TAL1 DP-like / TAL1 αβ-like |
| NKX2-5 | 7 | NA 6/7, ETP-like 1/7 | ETP-like, TAL1 αβ-like |
| NUP214 | 5 | MLLT10 3/5, HOXA9 TCR 1/5, NA 1/5 | MLLT10 |
| NUP98 | 6 | ETP-like;LMO2 γδ-like 2/6, MLLT10 2/6, KMT2A 1/6, NA 1/6 | ETP-like / MLLT10 |

### When you get `NA`
Possible reasons, roughly from most to least common:
1. **Rare subtype not in the model.** Check the top two vote fractions and the metagene plots. Confirm with fusion detection or genomic analysis (e.g. *NUP214*, *NUP98*, *NKX2-5* or *STAG2::LMO2* rearrangements).
2. **Heavy imputation.** Probabilities fall when many genes are imputed, especially for ETP-like and LMO2 γδ-like. If the app shows an imputation warning, re-quantify with a full annotation if you can.
3. **Low blast content or a non-T-ALL sample.** In 499 B-ALL samples, almost all calls were NA under the published thresholds. A sample dominated by normal cells can also fall below every threshold.
4. **Transitional or mixed state**, e.g. a relapse that changed subtype, or a sample with two leukaemic populations. Compare with the matched diagnosis sample if available.

### When `corrected` lists several subtypes (`A;B`)
More than one subtype passed its threshold. This happens mostly between transcriptionally adjacent subtypes (ETP-like / LMO2 γδ-like, NKX2-1 / TAL1 DP-like, MLLT10 / KMT2A / HOXA9 TCR). Report it as *ambiguous between A and B*, use the metagene plots to compare the two signatures, and resolve it with the genomic driver if known.

### ETP-like vs LMO2 γδ-like
These two subtypes share an immature program. TASC therefore uses a low threshold for LMO2 γδ-like (0.089), plus the override described in [Section 3](#3-understanding-the-output). An ETP-like call with a noticeable LMO2 γδ-like vote fraction is worth checking for *LMO2* rearrangements.

`[TO COMPLETE: any further rare-subtype tips, e.g. SPI1 fusions in TAL/LMO2-altered cases, KMT2A-rearranged B-ALL calling KMT2A/ETP-like]`

## 5. The TASC risk score

- **What it is:** a weighted combination of the 13 subtype vote fractions, after per-subtype recalibration. The weights come from a ridge-penalised Cox model of **overall survival** in the outcome-annotated patients of the reference cohort (Pölönen *et al.*, 2024). Subtypes that pushed the score up in the reference cohort include LMO2 γδ-like, TME-enriched and SPI1; MLLT10 pushes it down the most.
- **Risk classes:** Low, Intermediate and High follow the reference-cohort tertiles. They describe relative risk compared with that cohort, under its treatment protocol.
- **It is computed for every sample, including NA calls.** The score uses the whole probability profile rather than the single call. For rare subtypes, it ranked outcomes in the expected order in the reference cohort, but less strongly (Kendall τ = 0.33).
- **Use it together with MRD.** The score adds information to MRD rather than replacing it.

`[TO COMPLETE: citation of the TASC paper and external validation once published; contact for questions]`
