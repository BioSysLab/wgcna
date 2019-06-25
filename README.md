# wgcna
Code to perform WGCNA/Network based drug screening - Neidlin et al. 2019 

This repository includes all files to reproduce the work in Neidlin et al. 2019
More specifically the following analyses can be performed:

- Weighted gene co-expression network analysis : WGCNA folder 
=> Compute co-expressed and preserved modules across joint tissues. Extract driver genes of these modules

- Network based drug screening : NetworkAnalysis folder
=> Based on the driver genes (disease signature), compute drug-disease proximity measures in a background PPI network

- Gene set enrichment analysis : GSEA folder
=> Perform differential gene expression and gene set enrichment analysis with the piano package

-----Additional Scripts--------
- Data pre-processing with normalization and outlier removal : Preprocessing folder
- Stability analysis of the WGCNA algorithm : Stability folder
-------------------------------

Before running the R scripts, the following data needs to be downloaded:

Download from GEO
Synovium=GSE46750, Subchondral bone=GSE51588, Meniscus=GSE98918, Cartilage=GSE117999

Any questions/requests should be directed to:
Michael Neidlin

michael.neidlin@gmail.com

Athens, 25.06.2019
