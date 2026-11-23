## Project Overview

This project aims to computationally reproduce selected single-cell RNA sequencing (scRNA-seq) analyses presented in Figures 1 and 2 of the following research article:

*Integrative Single-Cell and Bulk Transcriptomic Analysis Identifies Macrophage-Related Gene Signatures Predictive of Hepatocellular Carcinoma in Cirrhosis.*

[Read the original publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC12562699/)

The original study's analytical methods and reported parameters were consulted as references to guide data preprocessing, integration, cell type annotation, and downstream macrophage analyses. This project represents an independent computational reproduction of selected analyses, rather than an exact reproduction of the entire published study.

## Datasets

Four publicly available scRNA-seq datasets from GEO were included:

| GEO Accession | Disease |
|:---|:---|
| [GSE136103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103) | Cirrhosis |
| [GSE186343](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186343) | HBV |
| [GSE242889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242889) | HCC |
| [GSE202642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202642) | HCC |

**Total number of cells analyzed:** Approximately 142,541.

## Analysis Workflow

### 1. Quality Control and Data Preprocessing

- Created Seurat objects from individual scRNA-seq datasets.
- Performed quality control (QC) to remove low-quality cells.
- QC thresholds were independently selected based on the characteristics of the datasets and informed by previously published studies.
  - `min.cells = 3`

  - `min.features = 300`

  - `max.mt = 25%`

  - `max.features` and `max.count` were individually determined for each dataset.

### 2. Data Integration and Batch Effect Correction

- Integrated the four datasets using the Reciprocal PCA (RPCA) integration approach implemented in Seurat.
- Performed reduction, clustering, and visualization to characterize cellular populations.

### 3. Differential Expression Analysis and Cell Type Annotation

- Differentially expressed genes (DEGs) across cell clusters were identified using the `FindAllMarkers()` function.

- Performed marker-based cell type annotation using a consensus approach by comparing cluster marker genes against three sources:

  - Canonical cell type markers curated from published literature.

  - PanglaoDB.

  - CellMarker 2.0.

  The annotation strategy was adapted from the following resource: [Single Cell RNA-seq — Marker-Based Cell Type Annotation](https://github.com/prpanigrahi/single_cell_annotation_course/tree/main/marker_based_annotation)

### 4. Macrophage Subset Analysis

- Extracted macrophage populations from the integrated dataset.
- Compared macrophage transcriptional profiles between HCC and cirrhosis groups.

#### Functional Enrichment Analysis

Conducted Gene Ontology (GO) enrichment analysis on differentially expressed genes between HCC and cirrhosis macrophages.

#### M1/M2 Macrophage Gene Signature Scoring

Compiled M1/M2 macrophage gene signatures from published literature and company references. Calculated macrophage polarization scores using Seurat's AddModuleScore() function to compare M1/M2 phenotypes across hepatitis, cirrhosis, and HCC.

## Repository Structure

``` text
project_aug/
│
├── README.md
│
├── code/
│   ├── packages.R
│   ├── functions.R
│   ├── 01_Preprocessing.qmd
│   ├── 02_Integration.qmd
│   ├── 03_Downstream_Analysis.qmd
│   ├── 04_Downstream_Analysis_2.qmd
│   └── 05_Downstream_Analysis_3.qmd
│
├── data/    #Download data here                     
│   ├── GSE136103_Data/
│   │   └── use/
│   ├── GSE186343_Data/
│   │   └── use/
│   ├── GSE202642_Data/
│   │   └── use/
│   └── GSE242889_Data/
│       └── use/
│
├── markers/
│   └── Cell_marker_Human.xlsx
│   └── PanglaoDB_markers_27_Mar_2020.tsv
│   └── Marker Summary.xlsx
|
├── results/
│   ├── pictures/
│   └── RDS/
│
└── session_info.txt
```

## References

1.  Zhang, Z.; Zeng, C.; Yong, X.; Zhou, W.; Xie, Y.; Shu, J. Integrative Single-Cell and Bulk Transcriptomic Analysis Identifies Macrophage-Related Gene Signatures Predictive of Hepatocellular Carcinoma in Cirrhosis . *Genes* **2025**, *16*, 1213. <https://doi.org/10.3390/genes16101213>

2.  Panigrahi PR. *Single Cell RNA-seq — Marker-Based Cell Type Annotation*. [GitHub repository](https://github.com/prpanigrahi/single_cell_annotation_course/tree/main/marker_based_annotation).

3.  Oscar Franzén, Li-Ming Gan, Johan L M Björkegren, PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data, *Database*, Volume 2019, 2019, baz046, <https://doi.org/10.1093/database/baz046>

4.  Congxue Hu, Tengyue Li, Yingqi Xu, Xinxin Zhang, Feng Li, Jing Bai, Jing Chen, Wenqi Jiang, Kaiyue Yang, Qi Ou, Xia Li, Peng Wang, Yunpeng Zhang, CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data, *Nucleic Acids Research*, Volume 51, Issue D1, 6 January 2023, Pages D870–D876, <https://doi.org/10.1093/nar/gkac947>
