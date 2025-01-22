## scRNA-Scanpy-WebApp  

### Overview  
This **Streamlit**-based web application provides an interactive platform for processing and analyzing single-cell RNA sequencing (scRNA-seq) data using **Scanpy**. Workflow includes data upload, quality control (QC), normalization, dimensionality reduction, clustering, and marker discovery.

---

### Features  
- **Data Upload & QC**  
  - Supports multiple formats: `10x Genomics (H5)`, `AnnData (H5AD)`, `CSV`, `MTX`, `Loom`, and Excel.  
  - Automated detection of mitochondrial, ribosomal, and hemoglobin genes for QC.  
  - Violin and scatter plots for visualizing quality metrics.

- **Processing**  
  - Doublet detection with Scrublet.  
  - Data normalization and log transformation.  
  - Highly variable gene selection (HVG).

- **Dimensionality Reduction**  
  - PCA for feature space reduction.  
  - Interactive sliders for parameter tuning.

- **Clustering & Visualization**  
  - UMAP and t-SNE integration.  
  - Leiden clustering with resolution adjustments.

- **Marker Discovery**  
  - Differential expression analysis for marker genes.  
  - Heatmaps and rank plots for visualizing key markers.

---
