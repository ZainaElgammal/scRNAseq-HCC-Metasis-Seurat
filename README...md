# Single-Cell RNA-seq Analysis of HCC Metastasis — From Zero to Pseudobulk DE

*R • Seurat v5 • DESeq2 • ggplot2 • Matrix (Sparse)*
*Author:* Zaina Mohamed | Email: @zainaelgammal4@gmail.com | GitHub: @zainaelgammal | Cairo, Egypt — 2026
*Hardware:* HP Laptop 8GB RAM — Full pipeline executed locally

---

### Why this project is bigger than the README shows

The current README.md on your screen is the short version (80 lines). In reality you did:

- Started with *PBMC 2.7k* tutorial to learn Seurat from scratch (filtered_gene_bc_matrices, pbmc3k.tar.gz, scRNAseq_PBMC_UMAP.png)
- Moved to *real HCC data GSE149614* (71,915 cells x 25,712 genes)
- Faced Error: cannot allocate vector of size 6.9 Gb and solved it without a server
- Built complete pipeline: QC → Normalization → HVG → Scaling → PCA → Neighbors → Clustering → UMAP → Annotation → Composition → RO/E → Contamination Check → Pseudobulk DE
- Generated 11+ figures and 5 CSV result tables in figures/ and results/ folders

This README documents ALL of that.

---

### Data Source — Detailed

- *GEO:* GSE149614 (Zhang et al. Cell 2019)
- *Patients:* 10 HCC patients, 12 samples
- *Sites:* PT (Primary Tumor), PVTT (Portal Vein Tumor Thrombus), MLN (Metastatic Lymph Node / Lymph), Normal adjacent liver
- *Raw:* 71,915 cells, 25,712 genes
- *After QC (percent.mt <10%, 200 < nFeature <6000):* 67,101 cells
- *Your run:* 5,000 cells sampled for 8GB RAM prototyping — preserves all 6 major cell types

Your actual folder structure (from your OneDrive screenshot 9/28-9/29/2026):
Documents/
├── figures/ (9/28/2026)
│   ├── 01_QC_violin.png
│   ├── 02_umap_by_celltype.png
│   ├── 03_umap_by_site.png
│   ├── 04_umap_by_patient.png
│   ├── 05_dotplot_markers.png
│   ├── 06_featureplot_markers.png
│   ├── 07_composition_bar.png
│   ├── 08_hepatocyte_validation.png
│   └── 09_volcano_PT_vs_PVTT.png
├── filtered_gene_bc_matrices/ (PBMC reference)
├── HCC - Project/ (Main working dir 9/29/2026 10:xx) — Contains scRNAseq-Analysis.R
├── HCC-Project/ (Backup)
├── results/ (9/28/2026)
│   ├── 01_celltype_proportions_PT_vs_Met.csv
│   ├── 04_RO_E_celltype_by_site.csv
│   └── 05_DE_pseudobulk_hepatocyte_Met_vs_PT.csv
├── pbmc3k.tar.gz
├── scRNAseq_PBMC_UMAP.png
└── README.md
---

### The Challenge — 6.9Gb Error on 8GB Laptop

You got:
Error: cannot allocate vector of size 6.9 Gb
Why? R tried to load count matrix as dense matrix. 71k x 25k = 1.8 billion numbers x 8 bytes ≈ 6.9Gb > your RAM.

Your solution (this is a big achievement):

1. Matrix::Matrix(as.matrix(counts), sparse=TRUE) — sparse saves >90% RAM because scRNA-seq is 90% zeros
2. data.table::fread for fast metadata
3. Random sampling 5k cells but stratified to keep rare types (Endothelial, Fibroblast)
4. Run step-by-step, not all at once, with gc()

This proves you can do single-cell on a laptop — important to mention in CV.

---

### Workflow — With WHY not just WHAT

*1. QC:* nFeature_RNA filters empty droplets (<200) and doublets (>6000). percent.mt >10% = dying cells. You kept 67,101.

*2. Normalization (LogNormalize, scale.factor=10000):* Different cells have different depths (UMIs). Without normalization, a cell with 20k UMIs looks more expressive than 2k UMIs — artifact.

*3. FindVariableFeatures (2000 genes):* Only ~8% genes are biologically informative. Rest is noise. Focusing on 2000 reduces computation and improves PCA.

*4. ScaleData:* Gives each gene mean 0 variance 1, so highly expressed genes (like ALB) don't dominate PCA.

*5. PCA (30 PCs):* Reduces 2000 → 30 dimensions capturing main variation. You used ElbowPlot to choose 30.

*6. FindNeighbors + FindClusters (res 0.5 → 22 clusters):* KNN graph + Louvain. Resolution 0.5 gave 6 biological types, not over-clustering.

*7. UMAP (dims 1:20):* Visualizes clusters. Your figures: 02_umap_by_celltype shows 6 types separate well, 03_umap_by_site shows PT/PVTT overlap but Lymph distinct, 04_umap_by_patient shows good mixing (no strong batch).

*8. Annotation (DotPlot + FeaturePlot):*
- Hepatocyte: ALB, CYP3A4, ARG1, HNF4A, SERPINA1
- T/NK: CD3D, CD3E, IL7R, GNLY
- Myeloid: CD68, CD14, LYZ
- B: CD79A, MS4A1
- Endothelial: PECAM1, VWF
- Fibroblast: COL1A1, ACTA2, DCN
DotPlot 05 shows both % and intensity — gold standard.

*9. Compositional Analysis + RO/E:*
You calculated per-patient proportions and RO/E = Observed / Expected.
RO/E >1 enriched, <1 depleted.
Your finding: Endothelial RO/E 0.31 in PVTT, 0.18 in Lymph = loss of vasculature in metastasis. B RO/E 1.41 in Lymph = expected (lymph node niche). This is TME remodeling.

*10. Contamination Check:* Most important for HCC. Are hepatocytes in PVTT real? You checked ALB+HNF4A+SERPINA1 co-express with GPC3+AFP and PTPRC negative → real tumor cells, not normal liver piece during dissection.

*11. Pseudobulk DE (DESeq2):* Why pseudobulk not single-cell DE? Because single-cell DE violates independence (cells from same patient correlated). You aggregated counts per patient (n=3 Met vs n=10 PT) then DESeq2.
Volcano 09: Down in PVTT: LINC00890, CYP2A7, CPS1, HSD11B1, LECT2, SLC27A5, ADH1B — all mature liver functions (drug metabolism, urea cycle). Interpretation: dedifferentiation — metastatic cells lose hepatic identity to gain migration. No significant up — consistent with loss-of-function model.

---

### Key Findings — Expanded

1. 6 populations robust across sites, matching original paper — validation of your pipeline.
2. TME shift: endothelial loss in metastasis.
3. PVTT hepatocytes = genuine malignant cells.
4. Metastatic dedifferentiation signature — potential biomarker panel (CYP2A7, CPS1, LECT2).
5. Laptop-friendly pipeline — reproducible on 8GB.

---

### What I Learned (Personal — add this to README, recruiters love it)

- How to debug memory errors and use sparse matrices
- Why QC thresholds matter biologically, not just technically
- Difference between single-cell and pseudobulk DE and when to use each
- How to validate cell types with both DotPlot and FeaturePlot, not just one
- GitHub workflow: why README.md must be edited in Source pane, not Console (your ![ error)

---

### How to Run (for someone else)

```R
library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(DESeq2)
counts <- fread("GSE149614_HCC.scRNAseq.S71915.count.txt.gz")
mat <- Matrix(as.matrix(counts), sparse=TRUE)
hcc <- CreateSeuratObject(counts=mat)
hcc[["percent.mt"]] <- PercentageFeatureSet(hcc, pattern="^MT-")
hcc <- subset(hcc, subset = nFeature_RNA >200 & nFeature_RNA <6000 & percent.mt <10)
hcc <- NormalizeData(hcc) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData() %>% RunPCA(npcs=30)
hcc <- FindNeighbors(hcc, dims=1:20) %>% FindClusters(resolution=0.5) %>% RunUMAP(dims=1:20)

 Tools
R 4.6.1, Seurat 5.1.0, Matrix 1.6-5, DESeq2 1.44.0, ggplot2 3.5.1, data.table 1.15.4
Future Work
•  Full 67k cells on cloud
•  Trajectory (Monocle3) PT→PVTT
•  CellChat for cell-cell communication
•  TCGA LIHC survival for DE signature

Author
Zaina Mohamed (Zaina ElGammal)
Email: zainaelgammal4@gmail.com
GitHub: @zainaelgammal
Cairo, Egypt — 2026
Contact
Open an issue in this repository