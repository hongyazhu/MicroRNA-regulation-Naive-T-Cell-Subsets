# miR-Inf
The miRNA target prediction framework, miR-Inf, is adapted from a transcriptional regulatory network, Inferelator: [Miraldi et al., Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells](https://genome.cshlp.org/content/early/2019/02/19/gr.238253.118). Code of this directory is adapted from the original code for Inferelator (https://github.com/emiraldi/infTRN_lassoStARS).

Key adaptaions made from the Inferelator algorithm to the miR-Inf framework:
1. Using estimated decay rates (see figure 2) instead of gene expression matrix as target gene levels. Estimated decay rates can better reflect miRNA regulation than gene expression levels.
2. Using miRNA expression levels (small RNA-seq) instead of TF expression levels as expression levels of regulators (in this case miRNAs). importGeneExpGeneLists_extraTFmRNA.m is adapted from importGeneExpGeneLists.m to import the miRNA expression matrix as NEW potential regulator levels.
3. Limiting solves to be positive only, as miRNAs can only accelerate decay of genes. When using function glmnetSet, set options.cl as instead of [-Inf;Inf] in getMLassoStARSlambdaRangePerGene.m and getMLassoStARSinstabilitiesPerGeneAndNet.m.
4. Considering only miRNAs with target sites in a gene as potential regulators. A miRNA is not considered for a certain gene if the value for that TF and gene is 0 in the input prior matrix, achieved by changes in estimateInstabilitiesTRNbStARS.m (descibed below).
5. Prior is not limited to 0 and 1. Strength of the prior is correlated with the strength of miRNA target prediction from Targetscan (which is reflected by total context score). If the absolute value of total context score is larger, the miRNA targeting is predicted to be stronger, and the prior will be more strongly incorporated. This is achieved by changes in estimateInstabilitiesTRNbStARS.m.


Here are the detailed changes in the scripts:

Adaptation 2:
adding infLassoStARS/importGeneExpGeneLists_extraTFmRNA.m 

Adaptation 3:
getMLassoStARSlambdaRangePerGene.m line 116 and getMLassoStARSinstabilitiesPerGeneAndNet.m line 74.
preivous: options = glmnetSet;
now: opts.cl = [0;Inf]; options = glmnetSet(opts);

Adaptation 4:
estimateInstabilitiesTRNbStARS.m
at row 135, add:
'
is_one = priorWeightsMat == 1;
priorWeightsMat( is_one ) = inf;'

Adaptation 5:
estimateInstabilitiesTRNbStARS.m row 114.
previous: priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*abs(sign(priorMat)); 
now: priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*priorMat;
Priors should be scaled to 0-1 before input.

