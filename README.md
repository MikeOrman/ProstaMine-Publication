# ProstaMine Shiny App

**ProstaMine** is a bioinformatics tool for identifying **subtype-specific co-alterations** associated with aggressiveness in prostate cancer.  
It integrates **genomic, transcriptomic, and clinical** data from multiple prostate cancer cohorts to prioritize gene co-alterations enriched in metastatic disease and associated with disease progression.

🧬 **Web app:** [https://bioinformatics.cuanschutz.edu/prostamine](https://bioinformatics.cuanschutz.edu/prostamine)  
📄 **Publication:** [Frontiers in Pharmacology (2024)](https://doi.org/10.3389/fphar.2024.1360352)

---

## 🚀 Overview

ProstaMine systematically mines prostate cancer data to:

- Identify **molecular subtype–specific co-alterations** (e.g., in *NKX3-1-loss* or *RB1-loss* tumors)  
- Integrate genomic, transcriptomic, and clinical outcomes  
- Rank alterations associated with **metastasis** and **biochemical relapse**  
- Visualize and export results via a **user-friendly R Shiny interface**

---

## Installation and Setup

You can run ProstaMine directly in your browser at the link above

## ⚙️ Using the App

### 🧩 Step 1 – Select a Molecular Subtype

Choose a genomic alteration to define your prostate cancer molecular subtype  
(e.g., **NKX3-1 loss**, **RB1 loss**, **PTEN loss**).  

ProstaMine automatically groups tumors into:

- **ST (Subtype):** samples with the selected alteration  
- **WT (Wild-type):** samples without the alteration  

---

### 🧠 Step 2 – Set Analysis Parameters

Adjust filtering parameters in the sidebar to control stringency.

| Parameter | Description | Default |
|------------|--------------|----------|
| **Primary co-alteration frequency difference** | Minimum frequency difference between ST vs WT primary tumors | 0.05 |
| **Metastatic co-alteration frequency difference** | Minimum frequency difference between ST metastatic vs ST primary tumors | 0.05 |
| **Differential gene expression FDR** | FDR threshold for transcriptomic hits | 0.20 |
| **Survival p-value cutoff** | p-value cutoff for progression-free survival association | 0.05 |

---

## 📊 Output and Visualization

### 🧬 1. Genomic Analysis
- Displays significant co-alterations enriched in the selected subtype  
- Presented as alteration heatmaps and frequency plots  

### 🧠 2. Transcriptomic Analysis
- Highlights genes with concordant differential expression in primary and metastatic tumors  
- Fold change and FDR values are reported  

### 🩸 3. Clinical Association
- Kaplan–Meier progression-free survival curves  
- Gleason grade enrichment tests  
- Identifies alterations linked to poor prognosis  

### 🧮 4. Ranked Co-alteration Table

All results are ranked using the **ProstaMine Score**:

```text
ProstaMine Score = 0.3 × (Primary ΔFreq Rank)
                 + 0.3 × (Metastatic ΔFreq Rank)
                 + 0.4 × (Survival Rank)
```

Higher ProstaMine scores indicate stronger association with aggressive disease.

## 💾 Exporting Results

Each results tab includes download options for:

- **Ranked co-alteration table** (`.csv`)
- **Summary plots** (`.png`, `.pdf`)
- **Parameter log** (for reproducibility)

---

## 📚 Data Sources

ProstaMine leverages harmonized multi-omic data via the  
[`curatedPCaData`](https://github.com/FIMM-CURATED/curatedPCaData) R package  
(*Laajala et al.,* *Scientific Data*, 2023), integrating genomic, transcriptomic, and clinical features from six key prostate cancer cohorts:

- **TCGA-PRAD**
- **Taylor et al., 2010**
- **Barbieri et al., 2012**
- **Baca et al., 2013**
- **Hieronymus et al., 2014**
- **Abida et al., 2019**

---

## 🧠 Citation

If you use **ProstaMine**, please cite:

> **Orman MV, Sreekanth V, Laajala TD, Cramer SD, and Costello JC.**  
> *ProstaMine: a bioinformatics tool for identifying subtype-specific co-alterations associated with aggressiveness in prostate cancer.*  
> **Frontiers in Pharmacology** (2024). doi:[10.3389/fphar.2024.1360352](https://doi.org/10.3389/fphar.2024.1360352)
