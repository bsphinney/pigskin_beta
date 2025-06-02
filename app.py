
# Save this as app.py and run with: streamlit run app.py

import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load full bidirectional dataset
full_data = pd.read_csv("All_Significant_Proteins_Log2FCgt0.58_or_lt_neg0.58_Qlt0.05.csv")
muscle_data = pd.read_csv("Muscle_Repair_Related_Proteins_Bidirectional.csv")

# Sidebar filters
st.sidebar.header("🔍 Filters")
show_muscle_only = st.sidebar.checkbox("Show only muscle-repair proteins", value=False)
direction = st.sidebar.selectbox("Regulation Direction", ["Upregulated", "Downregulated", "Both"])
log2fc_thresh = st.sidebar.slider("Log2FC Threshold (abs)", 0.58, 5.0, 0.58)
qval_thresh = st.sidebar.slider("Maximum Q-value", 0.0001, 0.05, 0.05)

# Apply filters
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

# Title
st.title("🧬 Wound Healing Proteomics Dashboard")
st.markdown("Explore significantly regulated proteins by direction, function, and treatment.")

# Heatmap
st.subheader("📊 Heatmap of Log2 Fold Changes")
pivot = filtered.pivot_table(index="Genes", columns="Comparison (group1/group2)", values="AVG Log2 Ratio")
if not pivot.empty:
    fig, ax = plt.subplots(figsize=(10, min(0.5 * len(pivot), 20)))
    sns.heatmap(pivot, cmap="coolwarm", center=0, annot=True, fmt=".2f", ax=ax, cbar_kws={"label": "Log2FC"})
    st.pyplot(fig)
else:
    st.info("No data to display with current filters.")

# Gene selector
if not filtered.empty:
    gene = st.selectbox("Select a gene to view details", sorted(filtered["Genes"].dropna().unique()))
    st.subheader(f"🔬 Details for {gene}")
    st.dataframe(filtered[filtered["Genes"] == gene])

# Show filtered data
st.subheader("🧾 Filtered Results")
st.dataframe(filtered)

# Download
st.download_button("Download Current Filtered Table", filtered.to_csv(index=False), file_name="filtered_proteins.csv")
if show_muscle_only:
    st.download_button("Download All Muscle Repair Proteins", muscle_data.to_csv(index=False), file_name="muscle_repair_proteins.csv")
