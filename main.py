import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import os
import tempfile
import gc
from io import BytesIO
import seaborn as sns

st.set_page_config(page_title="scRNA-Scanpy-WebApp", page_icon="🔬", layout="wide")


with st.sidebar:
    sidebar_container = st.container()
    
    with sidebar_container:
        st.header("🔬 Analysis Workflow")
        st.markdown("---")
        
        st.subheader("Main Steps:")
        st.markdown("""
        - **Data Upload & QC**  
        - **Processing & Normalization**  
        - **Dimensionality Reduction**  
        - **Clustering Analysis**  
        - **Marker Discovery**
        """)
        
        st.markdown("---")
        st.markdown("<br>", unsafe_allow_html=True)



# =====================================
# Core Functions 
# =====================================
def read_singlecell_data(file_path, format_hint=None):
    ext = os.path.splitext(file_path)[1].lower()
    
    try:
        if format_hint == '10x Genomics' or ext == '.h5':
            adata = sc.read_10x_h5(file_path)
            adata.var_names_make_unique()
            return adata
            
        elif format_hint == 'AnnData' or ext == '.h5ad':
            return sc.read_h5ad(file_path)
            
        elif ext == '.csv':
            return sc.read_csv(file_path)
            
        elif ext == '.mtx':
            adata = sc.read_mtx(file_path)
            adata.var_names_make_unique()
            return adata
            
        elif ext == '.loom':
            return sc.read_loom(file_path)
            
        elif ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path, engine='openpyxl')
            return sc.AnnData(df)
            
        elif ext in ['.txt', '.tsv']:
            delimiter = '\t' if ext == '.tsv' else None
            return sc.read_text(file_path, delimiter=delimiter)
            
        else:
            raise ValueError(f"Unsupported file format: {ext}")
            
    except Exception as e:
        raise ValueError(f"Error reading {file_path}: {str(e)}")

def plot_with_style(func, *args, **kwargs):
    sns.set_theme(style="whitegrid")
    fig = plt.figure(figsize=(8, 6), dpi=120)
    try:
        ax = func(*args, show=False, **kwargs)
        if ax is not None:
            plt.tight_layout()
            return ax.figure
        else:
            func(*args, show=False, ax=fig.gca(), **kwargs)
            plt.tight_layout()
            return fig
    except TypeError:
        func(*args, show=False, **kwargs)
        plt.tight_layout()
        return fig
    finally:
        plt.close(fig)


if 'adata' not in st.session_state:
    st.session_state.update({
        'adata': None,
        'processing_steps': [],
        'current_resolution': 0.5,
        'selected_markers': []
    })

# =====================================
# Main 
# =====================================
st.title("🔍 scRNA-Scanpy-WebApp ")
st.markdown("---")

with st.expander("📤 DATA UPLOAD", expanded=True):
    col1, col2 = st.columns([3, 1])
    with col1:
        uploaded_file = st.file_uploader(
            "Upload single-cell data", 
            type=['h5', 'h5ad', 'mtx', 'loom', 'csv', 'txt', 'xlsx'],
            help="Supports 10x Genomics, AnnData, CSV, Excel and more",
            label_visibility="visible"
        )
    with col2:
        file_format = st.selectbox(
            "File format", 
            ['Auto Detect', '10x Genomics', 'AnnData', 'CSV/TSV', 'Excel'],
            help="Manually specify if auto-detection fails"
        )
    
    st.markdown("""
    <div style="
        text-align: center;
        color: #666;
        margin-top: 1rem;
        font-size: 0.9em;
    ">
        Drag and drop files here<br>
        • Supported formats: H5, H5AD, MTX, LOOM, CSV, TXT, XLSX
    </div>
    """, unsafe_allow_html=True)
    
    if uploaded_file and st.button("🚀 Load Data", use_container_width=True):
        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = os.path.join(tmp_dir, uploaded_file.name)
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            
            try:
                with st.spinner('Analyzing data structure...'):
                    fmt = None if file_format == 'Auto Detect' else file_format
                    adata = read_singlecell_data(file_path, fmt)
                    
                    if 'sample' not in adata.obs:
                        sample_id = os.path.splitext(uploaded_file.name)[0]
                        adata.obs['sample'] = sample_id
                    
                    st.session_state.adata = adata
                    st.success(f"Loaded {adata.n_obs} cells with {adata.n_vars} features")
                    
            except Exception as e:
                st.error(f"Error loading file: {str(e)}")


if st.session_state.adata is not None:
    st.markdown("---")
    
    # QC
    with st.expander("🔍 QUALITY CONTROL", expanded=True):
        organism = st.selectbox("Select organism", ["Human", "Mouse"])
        mt_gene_prefix = "MT-" if organism == "Human" else "Mt-"
        
        with st.spinner('Calculating QC metrics...'):
            st.session_state.adata.var["mt"] = st.session_state.adata.var_names.str.startswith(mt_gene_prefix)
            st.session_state.adata.var["ribo"] = st.session_state.adata.var_names.str.startswith(("RPS", "RPL"))
            st.session_state.adata.var["hb"] = st.session_state.adata.var_names.str.contains("^HB[^(P)]")
            sc.pp.calculate_qc_metrics(
                st.session_state.adata, 
                qc_vars=["mt", "ribo", "hb"], 
                inplace=True, 
                log1p=True
            )

        viz_type = st.radio("Visualization type:", ["Violin Plots", "Scatter Plot"], horizontal=True)
        with st.spinner('Generating plot...'):
            if viz_type == "Violin Plots":
                fig, ax = plt.subplots(1, 3, figsize=(15, 5))
                for i, metric in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
                    sc.pl.violin(st.session_state.adata, metric, jitter=0.4, show=False, ax=ax[i])
                st.pyplot(fig)
                plt.close(fig)
            else:
                fig = plot_with_style(sc.pl.scatter, st.session_state.adata, 
                                    "total_counts", "n_genes_by_counts", color="pct_counts_mt")
                st.pyplot(fig)
                plt.close(fig)
            gc.collect()

        col1, col2 = st.columns(2)
        with col1:
            min_genes = st.slider("Minimum genes per cell", 0, 500, 100)
        with col2:
            min_cells = st.slider("Minimum cells per gene", 0, 10, 3)
        
        if st.button("🔧 Apply Filtering"):
            with st.spinner('Filtering...'):
                sc.pp.filter_cells(st.session_state.adata, min_genes=min_genes)
                sc.pp.filter_genes(st.session_state.adata, min_cells=min_cells)
            st.success(f"✅ Filtering applied! New shape: {st.session_state.adata.shape}")

    # processing 
    with st.expander("⚡ PROCESSING", expanded=False):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.checkbox("Detect Doublets"):
                with st.spinner('Running Scrublet...'):
                    if len(st.session_state.adata.obs['sample'].unique()) > 1:
                        sc.pp.scrublet(st.session_state.adata, batch_key="sample")
                    else:
                        sc.pp.scrublet(st.session_state.adata)
                    st.success("✅ Doublet detection complete!")
        
        with col2:
            if st.checkbox("Normalize Data"):
                with st.spinner('Normalizing...'):
                    st.session_state.adata.layers["counts"] = st.session_state.adata.X.copy()
                    sc.pp.normalize_total(st.session_state.adata)
                    sc.pp.log1p(st.session_state.adata)
                    st.success("✅ Normalization complete!")
        
        with col3:
            if st.checkbox("Select Features"):
                n_top_genes = st.number_input("HVG count", 500, 5000, 2000)
                with st.spinner('Selecting...'):
                    if len(st.session_state.adata.obs['sample'].unique()) > 1:
                        sc.pp.highly_variable_genes(
                            st.session_state.adata, 
                            n_top_genes=n_top_genes, 
                            batch_key="sample"
                        )
                    else:
                        sc.pp.highly_variable_genes(
                            st.session_state.adata, 
                            n_top_genes=n_top_genes
                        )
                    st.success("✅ Feature selection complete!")

    # dimensionality red.
    with st.expander("🌌 DIMENSIONALITY REDUCTION", expanded=False):
        if st.checkbox("Run PCA"):
            n_pcs = st.number_input("Number of PCs", 10, 100, 50)
            if st.button("Calculate PCA"):
                with st.spinner('Running PCA...'):
                    sc.tl.pca(st.session_state.adata, n_comps=n_pcs)
                    
                    plt.close('all')
                    sc.pl.pca_variance_ratio(
                        st.session_state.adata, 
                        n_pcs=n_pcs,
                        log=False,
                        show=False
                    )
                    
                    fig = plt.gcf()
                    fig.set_size_inches(8, 6)
                    fig.set_dpi(120)
                    plt.tight_layout()
                    st.pyplot(fig)
                    plt.close(fig)
                    gc.collect()

        st.subheader("UMAP Settings")
        n_neighbors = st.slider("Number of neighbors", 5, 100, 15)
        min_dist = st.slider("Minimum distance", 0.1, 1.0, 0.5)

        if st.checkbox("Run UMAP"):
            with st.spinner('Computing UMAP...'):
                sc.pp.neighbors(st.session_state.adata, n_neighbors=n_neighbors)
                sc.tl.umap(st.session_state.adata, min_dist=min_dist)
                fig = plot_with_style(sc.pl.umap, st.session_state.adata, color="sample")
                st.pyplot(fig)
                plt.close(fig)
                gc.collect()

    # clust.
    with st.expander("🧩 CLUSTERING & VISUALIZATION", expanded=True):
        st.subheader("Leiden Clustering")
        col1, col2 = st.columns([2, 1])
        
        with col1:
            resolution = st.slider(
                "Cluster resolution", 
                0.1, 2.0, 0.5,
                help="Higher values create more clusters"
            )
            if st.button("Run Leiden Clustering"):
                if 'neighbors' not in st.session_state.adata.uns:
                    st.warning("Compute neighbors first! Click 'Run UMAP' below")
                else:
                    with st.spinner('Clustering cells...'):
                        try:
                            sc.tl.leiden(
                                st.session_state.adata,
                                resolution=resolution,
                                key_added='leiden'
                            )
                            
                            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
                            
                            sc.pl.umap(
                                st.session_state.adata,
                                color='leiden',
                                ax=ax1,
                                show=False,
                                legend_loc='on data',
                                title=f'Leiden Clusters (res={resolution})',
                                frameon=False
                            )
                            
                            sc.pl.umap(
                                st.session_state.adata,
                                color='sample',
                                ax=ax2,
                                show=False,
                                title='Sample Distribution',
                                frameon=False
                            )
                            
                            plt.tight_layout()
                            st.pyplot(fig)
                            plt.close(fig)
                            st.success(f"Identified {len(st.session_state.adata.obs['leiden'].unique())} clusters")
                            
                        except Exception as e:
                            st.error(f"Clustering failed: {str(e)}")

        with col2:
            st.subheader("Cluster Sizes")
            if 'leiden' in st.session_state.adata.obs:
                cluster_counts = st.session_state.adata.obs['leiden'].value_counts().sort_index()
                st.dataframe(
                    cluster_counts.rename("Cell Count").to_frame(),
                    use_container_width=True,
                    height=300
                )

    # marker-analysis
    with st.expander("🔬 MARKER ANALYSIS", expanded=False):
        marker_genes = st.text_area(
            "Enter marker genes (JSON format)", 
            value='{"T Cells": ["CD3D", "CD8A"], "B Cells": ["CD19", "MS4A1"]}',
            height=150
        )
        
        if st.button("Analyze Markers"):
            try:
                marker_dict = json.loads(marker_genes)
                with st.spinner('Generating marker plot...'):
                    fig = plot_with_style(
                        sc.pl.dotplot,
                        st.session_state.adata,
                        marker_dict,
                        groupby='leiden',
                        standard_scale='var'
                    )
                    st.pyplot(fig)
                    plt.close(fig)
            except Exception as e:
                st.error(f"Error: {str(e)}")
        gc.collect()

    # differential-expression
    with st.expander("📈 DIFFERENTIAL EXPRESSION", expanded=False):
        if 'leiden' in st.session_state.adata.obs:
            de_cluster = st.selectbox(
                "Select cluster for DE analysis",
                st.session_state.adata.obs['leiden'].unique()
            )
            
            if st.button("Find DE Genes"):
                with st.spinner('Calculating differential expression...'):
                    sc.tl.rank_genes_groups(
                        st.session_state.adata, 
                        groupby='leiden', 
                        method='wilcoxon'
                    )
                    de_df = sc.get.rank_genes_groups_df(st.session_state.adata, group=de_cluster)
                    st.dataframe(de_df.head(20), use_container_width=True)
        gc.collect()


st.markdown("---")
with st.expander("📤 EXPORT RESULTS"):
    if st.session_state.adata is not None:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            st.session_state.adata.write(tmp.name)
            with open(tmp.name, "rb") as f:
                h5ad_data = f.read()
            
            st.download_button(
                "Download AnnData (.h5ad)",
                data=h5ad_data,
                file_name="analysis_results.h5ad",
                mime="application/octet-stream"
            )
        os.unlink(tmp.name) 

st.markdown("""
        ---
        **Created by Kristian Alikaj**  
        For more, visit [My GitHub](https://github.com/kris96tian) or [My Portfolio Website](https://kris96tian.github.io/)
""")
