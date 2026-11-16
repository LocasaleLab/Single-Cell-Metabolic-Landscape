# CCLE Colon (Large Intestine) RNA-seq: PCA, Clustering, and DESeq2

**Goal.** Explore biologically distinct groups within colon cancer cell lines using bulk RNA-seq, and identify gene programs that separate those groups.

---

## Dataset

- **Counts:** CCLE RNA-seq gene counts (`*.gct.gz`; genes × samples).
- **Metadata:** CCLE sample annotations (site, pathology, histology, sex, age, TCGA code, etc.).
- **Focus subset:** `Site_Primary == "large_intestine"` (58 samples).

---

## Methods (summary)

1. **Preprocessing**
   - Read GCT, keep columns as samples.
   - Filter to **protein-coding** genes via Ensembl/`biomaRt`.
   - Collapse **duplicate gene symbols** by summing counts; set `Description` as rownames; drop non-numeric cols.
   - Align counts ↔ metadata by sample intersection & order; clean metadata (`"missing-info"` for NA in categoricals).

2. **Exploration**
   - Build `DESeqDataSet`; filter low counts; **VST** transform.
   - **PCA** with `plotPCA()` and `prcomp()`; elbow plot; extract **gene loadings**.
   - **Hierarchical clustering** (top 150 most-variable genes) → 4 clusters observed; for DE kept **2 clusters** for balanced comparison.
   - **Heatmaps** with ComplexHeatmap + sample annotations.

3. **Differential Expression**
   - **DESeq2** using `design = ~ Cluster`; shrink LFC & visualize with **volcano plot**.

---

## Key findings (high-level)

- **PC1 (~40% var)** loads on **intestinal epithelial/differentiation markers** (e.g., *GPA33, CEACAM5, CDH17, LGALS4, PHGR1, PPP1R1B*).  
  → Right side: more differentiated/epithelial-like; left: de-differentiated.  
- **PC2 (~9% var)** highlights **proteolysis & metabolism/invasion** axes (e.g., *KLK6, KLK10, KRT23, SLC2A3, MSLN, ZBED2*).
- Unsupervised clustering using `DESeq2` split cell lines into two broad states:
  - **Cluster 1 (stromal/mesenchymal-leaning signals):** up in *SGIP1, PDPN, COL3A1, PTN, ANGPTL2, GDF6, KRTAP1-5, PODN, GNAT3, NTM*.
  - **Cluster 2 (epithelial/differentiated):** up in *SLC26A3, OLFM4, GPA33, FABP1, TINAG, UGT1A10, SPINK4, MUC17, PIGR, MOGAT2*.
- No single metadata field (sex, pathology, TCGA code) fully explained the variance—supporting a **transcriptional state** separation rather than a labeled covariate.

---

## Reproduce

**Requirements:** R (≥4.2), Quarto, and R/Bioc packages: `DESeq2`, `ComplexHeatmap`, `biomaRt`, `matrixStats`, `ggplot2`, `EnhancedVolcano`, `ggrepel`, `dplyr`, `tidyr`.

```bash
# from repo root
quarto render analysis/bulkRNAseq.qmd   # produces analysis/bulkRNAseq.pdf
```

If you don’t use Quarto, open the .qmd in RStudio and “Render”.

## Limitations & future directions

- Batch/confounders: Metadata incompleteness (e.g., age, race) limits covariate modeling; consider surrogate variable analysis (SVA) or RUV.

- Clustering k: Validate k with silhouette, gap statistic, and stability resampling; try consensus clustering.

- Pathway biology: Run GSEA/ORA on DEGs (Hallmarks/Reactome) to formalize epithelial vs mesenchymal/metabolic programs.

- Generalization: Cross-validate with TCGA COAD/READ and other colon cell-line panels.
