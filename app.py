
# Save this as app.py and run using: streamlit run app.py

import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from gprofiler import GProfiler
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform

# Load data
full_data = pd.read_csv("All_Significant_Proteins_Log2FCgt0.58_or_lt_neg0.58_Qlt0.05.csv")
muscle_data = pd.read_csv("Muscle_Repair_Related_Proteins_Bidirectional.csv")

# Sidebar filters
st.sidebar.header("🔍 Filters")
show_muscle_only = st.sidebar.checkbox("Show only muscle-repair proteins", value=False)
direction = st.sidebar.selectbox("Regulation Direction", ["Upregulated", "Downregulated", "Both"])
log2fc_thresh = st.sidebar.slider("Log2FC Threshold (abs)", 0.58, 5.0, 0.58)
qval_thresh = st.sidebar.slider("Maximum Q-value", 0.0001, 0.05, 0.05)

# Filter data
data = muscle_data if show_muscle_only else full_data

if direction == "Upregulated":
    filtered = data[(data["Direction"] == "Upregulated") & 
                    (data["AVG Log2 Ratio"] > log2fc_thresh) & 
                    (data["Qvalue"] < qval_thresh)]
elif direction == "Downregulated":
    filtered = data[(data["Direction"] == "Downregulated") & 
                    (data["AVG Log2 Ratio"] < -log2fc_thresh) & 
                    (data["Qvalue"] < qval_thresh)]
else:
    filtered = data[(abs(data["AVG Log2 Ratio"]) > log2fc_thresh) & 
                    (data["Qvalue"] < qval_thresh)]

# App title
st.title("🧬 Wound Healing Proteomics Dashboard")
st.markdown("Explore protein regulation across treatments with clustering, expression plots, and GO enrichment.")

# Clustered heatmap
st.subheader("📊 Clustered Heatmap of Log2 Fold Changes")
pivot = filtered.pivot_table(index="Genes", columns="Comparison (group1/group2)", values="AVG Log2 Ratio")

if not pivot.empty:
    # Hierarchical clustering for rows and columns using Euclidean distance and average linkage
    row_linkage = linkage(pdist(pivot.fillna(0), metric='euclidean'), method='average')
    col_linkage = linkage(pdist(pivot.fillna(0).T, metric='euclidean'), method='average')

    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)

    pivot_clustered = pivot.iloc[row_order, col_order]

    fig, ax = plt.subplots(figsize=(10, min(0.5 * len(pivot_clustered), 20)))
    sns.heatmap(pivot_clustered, cmap="coolwarm", center=0, annot=True, fmt=".2f", ax=ax, cbar_kws={"label": "Log2FC"})
    ax.set_title("Hierarchically Clustered Log2FC Heatmap")
    st.pyplot(fig)
else:
    st.info("No data to display with current filters.")

# Gene selector and expression plot
if not filtered.empty:
    gene = st.selectbox("Select a gene to view details", sorted(filtered["Genes"].dropna().unique()))
    gene_df = filtered[filtered["Genes"] == gene]

    st.subheader(f"🔬 Details for {gene}")
    st.dataframe(gene_df)

    st.subheader("📌 Expression Profile for Selected Gene Across Treatments")
    if not gene_df.empty:
        fig2, ax2 = plt.subplots()
        sns.barplot(data=gene_df, x="Comparison (group1/group2)", y="AVG Log2 Ratio", ax=ax2)
        ax2.set_title(f"Log2 Fold Change for {gene}")
        ax2.set_ylabel("Log2 Fold Change")
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha="right")
        st.pyplot(fig2)
    else:
        st.info("No data available for this gene.")

# GO enrichment
st.subheader("🧬 Gene Ontology Enrichment (GO:BP)")
go_genes = filtered["Genes"].dropna().unique().tolist()

if st.button("Run GO Enrichment"):
    gp = GProfiler(return_dataframe=True)
    results = gp.profile(organism='sscrofa', query=go_genes, sources=["GO:BP"])
    if results.empty:
        st.warning("No enrichment terms found.")
    else:
        top_results = results.sort_values("p_value").head(10)
        fig3, ax3 = plt.subplots(figsize=(8, 5))
        sns.barplot(data=top_results, y="name", x=-np.log10(top_results["p_value"]), ax=ax3)
        ax3.set_xlabel("-log10(P-value)")
        ax3.set_title("Top GO:BP Enrichment Terms")
        st.pyplot(fig3)
        st.dataframe(top_results)

# Download buttons
st.subheader("🧾 Filtered Results")
st.dataframe(filtered)
st.download_button("Download Filtered Results", filtered.to_csv(index=False), file_name="filtered_proteins.csv")
if show_muscle_only:
    st.download_button("Download All Muscle Repair Proteins", muscle_data.to_csv(index=False), file_name="muscle_repair_proteins.csv")
