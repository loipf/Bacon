# Bacon
Basic Atherosclerosis Comparative Online Network  
[a group project of bioinformatics students at the TUM]

Bacon provides a platform to analyze differential gene expression of microarray or RNA-Seq datasets, execute a gene set enrichment analysis and gives the possibility to compare interesting results to an already implemented atherosclerosis database with 13 (with different pipelines over 120) datasets with regards to atherosclerosis. Through ortholog comparison, it is possible to perform cross-species comparisons within the database which helps to evaluate the connection between human and mouse.

## pipeline how a user can analyse his data with Bacon
<div style="width:50%; height:50%">
![user example](./pic/use_case.png)
</div>


13 different datasets were analysed using following various mappers:
<div style="width:50%; height:50%">
![dataset pipeline](./pic/dge_pipeline.png)
</div>

---

## website with tools
[Rstudio server is necessary to host the apps]

![website home](./pic/website_home.png)

---
### differential gene expression analysis tool
![DGE prediction GUI](./pic/dge_prediction.png)
![DGE analysis GUI](./pic/dge_analysis.png)
![DGE qc GUI](./pic/dge_qc_plots.png)


---
### gene set enrichment analysis against 13 implemented datasets
![GSE database GUI](./pic/gse_db_overlap.png)

---
### ortholog comparison of mouse - human atherosclerosis genes
![ortholog comparison GUI](./pic/ortholog_comparison.png)


---
### the internal database can be extended by directly loading new datasets from GEO:
![GEO load GUI](./pic/load_GEO_data.png)





