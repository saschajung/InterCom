# InterCom - Reconstruction of functional cell-cell communication networks 

Reconstruct functional cell-cell communication networks from single-cell RNA-seq data of a tissue.

## Required software
  - R v3.6 or greater
  - textshape v1.6.0
  - ggplot2 v 3.3.0
  - stuRpkg v1.3
  - dplyr v0.8.5
  - doParallel v1.0.15
  - stringr v1.4.0
  - plyr 1.8.6
  - igraph v1.2.5
  - Matrix v1.2-18
  - reshape2 v1.4.4
  - RSpectra v0.16-0
  - snow v0.4-3
  - taRifx v1.0.6.2
  - gtools v3.8.2
  - data.table v1.12.8
  - rlist v0.4.6.1

## Running InterCom
The complete InterCom workflow can be invoked through a single command after sourcing the file "intercom_main_coexp.R"
```R
source("intercom_main_coexp.R")

InterCom(data,
         anno.tbl,
         species,
         sighot.cutoff=0.1,
         sighot.percentile=70,
         consv.thrs=0.05,
         ncores=4,
         sig.cutoff=0.9,
         z.score.cutoff=2,
         tissue.name,
         temp.folder.name,
         out.path
        )
```

## InterCom parameters
The main InterCom function comes with a variety of parameters for fine-tuning the output. However, in any case, the standard parameters work well in most of the cases. Below we describe every parameter, their ranges and standard values.

  - "data": scRNA-seq data of a tissue with celltype information as column names and row names in form of gene symbols. The input can be (normalized) count or TPM data.
  - "anno.tbl": Annotation data frame for the data with cell ids in first column and cell types in second. Importantly, the cell id column has to be named "cell.name" and the cell type information "cell.type". 
  - "species": Species information (HUMAN or MOUSE)
  - "sighot.cutoff": Cutoff parameter for the activation probability of high-probability intermediates in SigHotSpotter. Can be in the range of 0 to 1. (default: 0.1).
  - "sighot.percentile": Percentile parameter for SigHotSpotter. Can be in the range of 0 to 100 (default: 70).
  - "consv.thrs": Conservation threshold of TFs, ligands and receptors. Can be in the range of 0 to 1 (default: 0.05).
  - "ncores": Number of cores used for parallel execution. Can be any integer value between 1 to the number of available cores (default: 4).
  - "sig.cutoff": Significance cutoff for interactions. Can be between 0 (weakest) and 1 (strictest) (default: 0.9).
  - "z.score.cutoff": Cutoff parameter to determine significant associations between receptors and interface TFs. The value represents the number of standard deviations above (if greater 0) or below (if smaller than 0) the mean. Can be any numeric value. 
  - "tissue.name": A name for the run. Can be any string.
  - "temp.folder.name": A temporary folder. Can be any directory.
  - "out.path": The output folder. Can be any directory. If it does not exist, it will be created.

## InterCom output
The main output of InterCom is an RData file called "output_<tissue.name>.RData". The final interactome is stored in the object "output" as a data frame in the field "final".


