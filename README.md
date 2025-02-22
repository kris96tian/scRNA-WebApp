## scRNA-Scanpy-WebApp  

Visit App:  https://scrna-webapp.streamlit.app/

<img src="blob:chrome-untrusted://media-app/d7244cbd-044a-4cb8-b393-27546b9d6773"/>![image](https://github.com/user-attachments/assets/4fa04a1b-8be2-4d9d-b394-9565e02cf60d)


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
