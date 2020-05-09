# Welcome
This repository contains analysis scripts for Abbasi, Sinha, Labit et al. 2020 (Under Review). Single-cell datasets can be visualized on our [Wound Atlas](http://www.biernaskielab.ca/wound_atlas/).
![](images/2.%20Graphical%20Abstract.jpeg)

# Summary

## Abstract
Dermal fibroblasts exhibit considerable heterogeneity during homeostasis and in response to injury. Defining the lineage origins of wound healing fibroblasts and the regulatory programs that drive fibrosis or conversely promote regeneration will be essential toward developing therapeutics to improve wound healing outcomes. Using complementary fate mapping approaches, we show that hair follicle mesenchymal progenitors exhibit limited contribution to skin wound repair.  In contrast, extrafollicular dermal progenitors and their progeny marked by the quiescence-associated factor Hic1 generated the bulk of reparative fibroblasts and exhibited functional divergence corresponding to their location within wound neodermis; regeneration centrally and scar-forming in the periphery. Single-cell RNA-Seq of wound-responsive fibroblasts revealed unique transcriptional signatures (Crabp1, Fabp5, Prss35), regulatory network activity (Retinoic Acid Receptors, Runx1, Hox transcription factors), and epithelial-mesenchymal interactions that enabled competence for regeneration. Integration with scATAC-Seq further highlighted specific changes in chromatin accessibility within regeneration-associated loci. Finally, pharmacological modulation of RUNX1 and Retinoic Acid signaling or genetic deletion of Hic1 within wound-activated fibroblasts was sufficient to modulate healing outcomes, suggesting that reparative dermal fibroblast exhibit latent, but modifiable, regenerative capacity.

## Highlights:
1. Hair follicle dermal stem cells (hfDSCs) exhibit limited contribution to wound healing and hair follicle (HF) neogenesis.
2. Extrafollicular Hic1-lineage progenitors regenerate the injured dermis and populate neogenic HFs.
3. Distinct transcriptional and epigenetic changes enable fibroblasts to adopt divergent fates post-injury
4. Runx1, Retinoic Acid, and Hic1 regulate mesenchymal regenerative competence

# Data

## Single-cell RNA-Seq
NCBI GEO: [GSE108677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108677) <br/>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108677/suppl/GSE108677_RAW.tar
tar -xvf GSE108677_RAW.tar
```
NCBI SRA: SRP130923 <br/>
```
source activate sratoolkit
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
## Specify SRR_ID - obtained using SRA Run selector.
```

## Single-cell ATAC-Seq
NCBI GEO: [GSE131600](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131600) <br/>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131600/suppl/GSE131600_RAW.tar
tar -xvf GSE131600_RAW.tar
```
NCBI SRA: SRP199132 <br/>
```
source activate sratoolkit
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
## Specify SRR_ID - obtained using SRA Run selector.
```

## Single-cell MP datasets used for cross-tissue integration
Muscle NCBI GEO: [GSM2976778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2976778) (described in [Scott et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/31809738))<br/>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2976nnn/GSM2976778/suppl/GSM2976778_qsnt_barcodes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2976nnn/GSM2976778/suppl/GSM2976778_qsnt_genes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2976nnn/GSM2976778/suppl/GSM2976778_qsnt_matrix.mtx.gz
```
Heart NCBI GEO: [GSM2976778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2976778) (described in [Soliman et al. 2020](https://www.ncbi.nlm.nih.gov/pubmed/31978365))<br/>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4216nnn/GSM4216418/suppl/GSM4216418_Hic1tdTomato_undamaged_barcodes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4216nnn/GSM4216418/suppl/GSM4216418_Hic1tdTomato_undamaged_genes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4216nnn/GSM4216418/suppl/GSM4216418_Hic1tdTomato_undamaged_matrix.mtx.gz
```

# Toolkits used
`Seurat v.2.3.0` - Gene expression analysis shown in Figure 4 (`UpdateSeuratObject` used for minor re-analysis in Seurat v3). <br/>
`Seurat v.3.0.0` - scATAC-Seq + scRNA-Seq shown in Figure 5. <br/>
`Seurat v.3.0.0` - Cross-tissue integration shown in Figure 2. <br/>

# Contact
Dr. Jeff Biernaskie (jabierna@ucalgary.ca)<br/>
Sarthak Sinha (sarthak.chinoo@gmail.com)
