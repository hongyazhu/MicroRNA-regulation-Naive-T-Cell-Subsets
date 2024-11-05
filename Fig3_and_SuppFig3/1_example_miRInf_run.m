%% Example miR-Inf run
%% miR-Inf is adapted from Inferelator: Miraldi et al. (2018) "Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells"
%% code was adapted from example_workflow_Th17 by Dr. Emily Miraldi (https://github.com/emiraldi/infTRN_lassoStARS)

function 1_example_miRInf_run()
    tfaOpts = {'', '_TFmRNA'};
    for tfaOptNum = length(tfaOpts) % Only using '_TFmRNA' option
        for meanEdgesPerGene = [2] % reaches the largest network
            for seed = [7,18,99,26,57] % runs five times with five different seeds
                for lambdaBias = [0.05] % strong reinforcement
                    rng(seed)
                    tfaOpt = tfaOpts{tfaOptNum};
                    name_this_run = ['bias' num2str(100*lambdaBias) tfaOpt '_modelSize' num2str(meanEdgesPerGene)];
                    diaryFile = ['outputs/mirnaSeed_tsonly_noSP8_wt_ct02_seed' num2str(seed) '/' name_this_run '.log'];
                    matlabDir = 'miR-Inf';
                    outputFolder = ['outputs/mirnaSeed_tsonly_noSP8_wt_ct02_seed' num2str(seed) '/' name_this_run ];
                    normGeneExprFile = 'inputs/exonIntron_matrix/decay_rates_match_noNA_noInf_noall0_10exon.txt';
                    targGeneFile = 'inputs/exonIntron_matrix/targGene_match_noNA_noInf_noall0_10exon_tsmirnasexpressed1000_ct02.txt';
                    potRegFile = 'inputs/mirna_DE/mirnas_expressed1000.txt';
                    NEWpotRegmRNAlevels = 'inputs/mirna_exp_matrix/counts_vst_blindT_combat_expressed1000.txt';
                    tfaGeneFile = '';
                    priorName = 'ts_ct02_contextscore0to1_targGeneonly';
                    priorFile = 'inputs/prior/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly.tsv';
                    priorMergedTfsFile = '';
                    mkdir(outputFolder)
                    workflow_withArguments(diaryFile, matlabDir, outputFolder, ...
                        normGeneExprFile, targGeneFile, potRegFile, NEWpotRegmRNAlevels, tfaGeneFile, ...
                        lambdaBias, tfaOpt, priorName, priorFile, priorMergedTfsFile, ...
                        meanEdgesPerGene)
                 end
            end
        end
    end
end

% arguments - 
% matlabDir: directory of auxiliary matlab functions (infLassoStARS, glmnet, customMatlabFxns) [..]

% outputFolder: output directory. [./outputs]
% normGeneExprFile: normalized gene expression profiles [./inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt]
% targGeneFile: target genes names [./inputs/targRegLists/targetGenes_names.txt]
% potRegFile: potential regulators [./inputs/targRegLists/potRegs_names.txt]
% tfaGeneFile: genes for TFA [./inputs/targRegLists/genesForTFA.txt]

% tfaOpt: TFA estimated by equation or mRNA [ '_TFmRNA' or '' ]
% lambdaBias: prior inforcement strength [lambdaBiases = [1 .5 .25]; % correspond to no, moderate, and strong prior reinforcement]

% priorName: ['ATAC_Th17']
% priorFile: ['./inputs/priors/' priorName '.tsv']; 
% priorMergedTfsFile: ['./inputs/priors/' priorName '_mergedTfs.txt']

% meanEdgesPerGene: 15



function [] = workflow_withArguments(diaryFile, matlabDir, outputFolder, normGeneExprFile, targGeneFile, potRegFile, NEWpotRegmRNAlevels, tfaGeneFile, lambdaBias, tfaOpt, priorName, priorFile, priorMergedTfsFile, meanEdgesPerGene, gsFile, prNickName, prTargGeneFile)

diary(diaryFile)
diary on

maxNumCompThreads(4) 

%% example_workflow_Th17
% Use mLASSO-StARS to build a TRN from gene expression and prior
% information in four steps. Please refer to each function's help
% annotations for descriptions of inputs, outputs and other information.
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
% (2) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
% (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: March 29, 2018

%clear all %%% commented the two commands!
close all
%restoredefaultpath

%argument matlabDir = '..';

addpath(fullfile(matlabDir,'infLassoStARS'))
addpath(fullfile(matlabDir,'glmnet'))
addpath(fullfile(matlabDir,'customMatlabFxns'))

%% 1. Import gene expression data, list of regulators, list of target genes
% into a Matlab .mat object
%argument geneExprTFAdir = './outputs/processedGeneExpTFA';
geneExprTFAdir = append(outputFolder, '/processedGeneExpTFA'); % added!!!

mkdir(geneExprTFAdir)
%argument normGeneExprFile = './inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt';
%argument targGeneFile = './inputs/targRegLists/targetGenes_names.txt';
%argument potRegFile = './inputs/targRegLists/potRegs_names.txt';
%tfaGeneFile = 'inputs/gene_exp_matrix/counts_vst_blindT_combat.txt';
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

disp('1. importGeneExpGeneLists.m')
importGeneExpGeneLists_extraTFmRNA(normGeneExprFile,targGeneFile,potRegFile,NEWpotRegmRNAlevels,...
    tfaGeneFile,geneExprMat)

%% 2. Given a prior of TF-gene interactions, estimate transcription factor 
% activities (TFAs) using prior-based TFA and TF mRNA levels
%argument priorName = 'ATAC_Th17';
%argument priorFile = ['./inputs/priors/' priorName '.tsv']; % Th17 ATAC-seq prior
edgeSS = 50;
minTargets = 3;
[xx, priorName, ext] = fileparts(priorFile);
tfaMat = fullfile(geneExprTFAdir,[priorName '_ss' num2str(edgeSS) '.mat']);

disp('2. integratePrior_estTFA.m')
integratePrior_estTFA(geneExprMat,priorFile,edgeSS,...
     minTargets, tfaMat)

%% 3. Calculate network instabilities using bStARS

%argument lambdaBias = .5;
%argument tfaOpt = ''; % options are '_TFmRNA' or ''
totSS = 50;
targetInstability = .05;
lambdaMin = .01;
lambdaMax = 1;
extensionLimit = 1;
totLogLambdaSteps = 25; % will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5;
subsampleFrac = .63;
leaveOutSampleList = '';
leaveOutInf = '';
instabilitiesDir = fullfile(outputFolder,strrep(['instabilities_targ' ... %%% changed to outputFolder here!!!
    num2str(targetInstability) '_SS' num2str(totSS) leaveOutInf '_bS' num2str(bStarsTotSS)],'.','p'));
mkdir(instabilitiesDir)
netSummary = [priorName '_bias' strrep(num2str(100*lambdaBias),'.','p') tfaOpt];
instabOutMat = fullfile(instabilitiesDir,netSummary);

disp('3. estimateInstabilitiesTRNbStARS.m')
estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,...
    subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)

%% 4. For a given instability cutoff and model size, rank TF-gene
% interactions, calculate stabilities and network file for jp_gene_viz
% visualizations
%argument priorMergedTfsFile = ['./inputs/priors/' priorName '_mergedTfs.txt'];
try % not all priors have merged TFs and merged TF files
    ls(priorMergedTfsFile) 
catch
    priorMergedTfsFile = '';
end
%argument meanEdgesPerGene = 15;
targetInstability = .05;
networkDir = strrep(instabilitiesDir,'instabilities','networks');
instabSource = 'Network';
mkdir(networkDir);
networkSubDir = fullfile(networkDir,[instabSource ...
    strrep(num2str(targetInstability),'.','p') '_' ...
    num2str(meanEdgesPerGene) 'tfsPerGene']);
mkdir(networkSubDir)
trnOutMat = fullfile(networkSubDir,netSummary);
outNetFileSparse = fullfile(networkSubDir,[netSummary '_sp.tsv']);
networkHistDir = fullfile(networkSubDir,'Histograms');
mkdir(networkHistDir)
subsampHistPdf = fullfile(networkHistDir,[netSummary '_ssHist']);

disp('4. buildTRNs_mLassoStARS.m')
buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile,...
    meanEdgesPerGene,targetInstability,instabSource,subsampHistPdf,trnOutMat,...
    outNetFileSparse)


diary off

end
