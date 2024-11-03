# miR-Inf
Our miRNA target prediction framework, miR-Inf, is adapted from a transcriptional regulatory network, Inferelator (Miraldi et al., 2019; https://github.com/emiraldi/infTRN_lassoStARS).
Key adaptaions made:
* 1. Using estimated decay rates (see figure 2) instead of gene expression matrix as target gene levels. Estimated decay rates can better reflect miRNA regulation than gene expression levels.
* 2. Using miRNA expression levels (small RNA-seq) instead of TF expression levels as expression levels of regulators (in our case miRNAs). importGeneExpGeneLists_extraTFmRNA.m is adapted from importGeneExpGeneLists.m to import the miRNA expression matrix as NEW potential regulator levels.
* 3. Limiting solves to be positive only, as miRNAs can only accelerate decay of genes. When using function glmnetSet, set options.cl as instead of [-Inf;Inf] in getMLassoStARSlambdaRangePerGene.m and getMLassoStARSinstabilitiesPerGeneAndNet.m.
* 4. Considering only miRNAs with target sites in a gene as potential regulators. A miRNA is not considered for a certain gene if the value for that TF and gene is 0 in the input prior matrix, achieved by changes in estimateInstabilitiesTRNbStARS.m (descibed below).
* 5. Prior is not limited to 0 and 1. Strength of the prior is correlated with the strength of miRNA target prediction from Targetscan (which is reflected by total context score). If the absolute value of total context score is larger, the miRNA targeting is predicted to be stronger, and the prior will be more strongly incorporated. This is achieved by changes in estimateInstabilitiesTRNbStARS.m.


Here are the detailed changes in the scripts:

Adaptation 2:
infLassoStARS/importGeneExpGeneLists_extraTFmRNA.m 

Adaptation 3:
getMLassoStARSlambdaRangePerGene.m line 116 and getMLassoStARSinstabilitiesPerGeneAndNet.m line 74
preivous: options = glmnetSet;
now: opts.cl = [0;Inf]; options = glmnetSet(opts);

Adaptation 4:
estimateInstabilitiesTRNbStARS.m
at row 135, add:
'
is_one = priorWeightsMat == 1;
priorWeightsMat( is_one ) = inf;'

Adaptation 5:
estimateInstabilitiesTRNbStARS.m row 114 
previous: priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*abs(sign(priorMat)); 
now: priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*priorMat;
Priors should be scaled to 0-1 before input.


# README content for the Inferelator framework

# infTRN_lassoStARS

This repository contains a workflow for inference of transcriptional regulatory networks (TRNs) from gene expression data and prior information, as described in:

[Miraldi et al., Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells](https://genome.cshlp.org/content/early/2019/02/19/gr.238253.118).

From gene expression data and tables of prior information, the [example Th17 workflow](Th17_example/example_workflow_Th17.m) can be used to infer a TRN using modified LASSO-StARS, and relies upon [GlmNet in MATLAB](https://web.stanford.edu/~hastie/glmnet_matlab/index.html) to solve the LASSO. Workflow also includes TRN model evaluation based on precision-recall and ROC.

The resulting network can be visualized with TRN visualization software: [jp_gene_viz](https://github.com/simonsfoundation/jp_gene_viz).

Additional workflows are provided for:
* [Construct a prior transcriptional regulatory network from ATAC-seq data](priorConstruction/readme.md)
* [TF-TF module analysis: Discovery of TFs that co-regulate gene pathways](Th17_example/example_Th17_tfTfModules.m)
* [Identify "core" TF regulators for a subset of conditions or celltypes in the gene expression dataset](scTRN)
* [Gene-set enrichment analysis (GSEA) of a TF's positive and repressed target genes](scTRN)
* [Visualize TF degree, positive and negative edges, and overlap with TF-gene interactions in the prior](scTRN/viz_TF_degree.m)
* [Out-of-sample gene expression prediction](Th17_example/example_workflow_Th17_r2Pred.m), including calculation of R<sup>2</sup><sub>pred
* [Modeling of time-series gene expression with linear differential equations](Th17_example/example_workflow_Th17_timeLag.m), as in [Bonneau et al. (2006) Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2006-7-5-r36)

NOTE: For Mac users of MATLAB 2016b or later versions, you might need to install gfortran. We recommend the following: Install [Homebrew](https://brew.sh/), and then install gfortran with the Terminal command: ```"brew cask install gfortran"```.