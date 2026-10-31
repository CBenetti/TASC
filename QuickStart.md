# How to run TASC

## Using the app

Visit [TASC](https://01971d1c-d81d-a598-4238-5afb7b3e381a.share.connect.posit.cloud/) deployed version and immediately obtain your predictions.

### What you'll need:
A count matrix file, where transcriptomic counts from your samples of interest are stored.
Results from the most common alignment tools ([STAR](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/STAR_RNAseq.htm), [Salmon](https://combine-lab.github.io/salmon/)) are supported, as long as the counts are aggregated per gene, and genes are in [Ensembl ID](https://www.ebi.ac.uk/training/online/courses/ensembl-browsing-genomes/navigating-ensembl/investigating-a-gene/) format.

Alternatively, the general go-to format is a .tsv table with:
- Header line containing a "GeneID" field, followed by sample names
- Count matrix, with Ensembl IDs in the first column.

_Example_ :
| GeneID           | Sample_A | Sample_B | Sample_C | Sample_D | Sample_E |
|------------------|----------|----------|----------|----------|----------|
| ENSG00000139618  | 523      | 610      | 489      | 712      | 634      |
| ENSG00000157764  | 342      | 298      | 410      | 365      | 389      |
| ENSG00000141510  | 275      | 320      | 305      | 290      | 310      |
| ENSG00000155657  | 198      | 210      | 225      | 205      | 215      |
| ENSG00000171862  | 150      | 160      | 155      | 165      | 158      |

Missing genes required for the model are automatically imputated from an internal database based on [Polonen et.al](https://doi.org/10.1038/s41586-024-07807-0) pubblication.  
However, if more than 10% are missing, predictions could be unreliable.

In feature plots of metagenes (differentially exressed genes for each subtype) the required genes _will not be imputated_.
## Run the app on your machine
### Step 1
Clone the repo folder:

```
git clone --depth 1 https://github.com/CBenetti/TASC.git ./TASC
```

### Step 2

Opern an R or RStudio session, and make sure you have [Shiny](https://shiny.posit.co/) package installed.

### Step 3
In your session, run the app:

```
shiny::runApp("path_to_TASC_repo/TASC.R")
```
