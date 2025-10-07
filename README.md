ProstaMine Shiny App – README
Overview
ProstaMine is a bioinformatics tool for identifying subtype-specific co-alterations associated with aggressive features in prostate cancer.
It integrates genomic, transcriptomic, and clinical data across multiple prostate cancer cohorts and enables researchers to identify and rank co-altered genes associated with metastasis and biochemical relapse.
This README provides instructions for running the ProstaMine Shiny web application at: https://bioinformatics.cuanschutz.edu/prostamine
 
1. Getting Started
Accessing the App
No installation is required. Visit the Shiny app URL directly in your browser.
If running locally:
# Install required packages
install.packages(c("shiny", "tidyverse", "ComplexHeatmap", "ConsensusClusterPlus", 
                   "survival", "survminer", "karyoploteR"))
# Launch locally from GitHub
shiny::runGitHub("MikeOrman/ProstaMine-Publication", subdir = "ShinyApp")
 
2. Input Selection
Step 1: Choose Molecular Subtype
•	Select a molecular alteration (e.g., NKX3-1 loss, RB1 loss, PTEN loss) to define your subtype of interest.
•	ProstaMine automatically groups samples into:
o	ST (Subtype): tumors with the selected alteration
o	WT (Wild-type): tumors without the alteration
 
3. Analysis Parameters
Use the sidebar to set filtering parameters that control stringency:
Parameter	Description	Default
Primary co-alteration frequency difference	Minimum frequency difference between ST and WT primary tumors	0.05
Metastatic co-alteration frequency difference	Minimum frequency difference between ST metastatic and ST primary tumors	0.05
Differential gene expression FDR	FDR threshold for transcriptomic hits	0.2
Survival p-value cutoff	p-value cutoff for progression-free survival association	0.05
 
4. Output Tabs
1️ Genomic Analysis
•	Displays significant co-alterations enriched in the selected subtype.
•	Visualized as a heatmap and frequency plot.
2️ Transcriptomic Analysis
•	Shows genes with concordant differential expression in primary and metastatic tumors.
•	Fold change and FDR are displayed.
3️ Clinical Association
•	Kaplan–Meier survival plots and Gleason grade enrichment tests.
•	Identifies alterations linked to poor prognosis.
4️ Ranked Co-alteration Table
•	All hits are scored using the ProstaMine score:
•	ProstaMine Score = 0.3*(Primary ΔFreq Rank) + 0.3*(Metastatic ΔFreq Rank) + 0.4*(Survival Rank)
•	Higher scores indicate stronger association with aggressiveness.
 
5. Exporting Results
Each results tab includes download buttons to export:
•	Ranked co-alteration tables (.csv)
•	Summary plots (.png or .pdf)
•	Log of parameter settings used for reproducibility.
 
6. Data Sources
All underlying data are obtained from the curatedPCaData R package (Laajala et al., Sci Data, 2023), integrating genomic, transcriptomic, and clinical features from:
•	TCGA-PRAD
•	Taylor et al., 2010
•	Barbieri et al., 2012
•	Baca et al., 2013
•	Hieronymus et al., 2014
•	Abida et al., 2019
 
7. Citation
If you use ProstaMine, please cite:
Orman MV, Sreekanth V, Laajala TD, Cramer SD, and Costello JC (2024).
ProstaMine: a bioinformatics tool for identifying subtype-specific co-alterations associated with aggressiveness in prostate cancer.
Front Pharmacol. 15:1360352. doi: 10.3389/fphar.2024.1360352

<img width="468" height="661" alt="image" src="https://github.com/user-attachments/assets/b62ef724-7be3-4d31-971d-f6bb3ebb92e5" />
