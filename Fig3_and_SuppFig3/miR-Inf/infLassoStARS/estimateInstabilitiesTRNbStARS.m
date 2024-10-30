function estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,subsampleFrac,...
    instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit) 
% edit 1 is the addition of '
% is_one = priorWeightsMat == 1;
% priorWeightsMat( is_one ) = inf;'
% at row 135 (currently row 143 with my new comments) in estimateInstabilitiesTRNbStARS.m by Dr. Emily Miraldi.
% goal of this change: when predicting a certain gene, excluding all miRNAs that do not have target site at this gene (predicted by TargetScan) 
% (inf values will be excluded as shown at row 39-41 in getMLassoStARSlambdaRangePerGene.m)
% edit 2 is modification on row 114 (currently row 122 with my new comments)
% goal is to make the strength of prior correlate with TargetScan total context score (not just 0 or 1)
%% estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
%     totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,subsampleFrac,...
%     instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)
%% estimateInstabilitiesMLassoStars(path2data,priorName,priorTfaFile,lambdaBiases,...
%     tfaOpt,cvSubFolder,totSS,kfoldCvs)
%% Goal: Estimate mLASSO-StARS instabilities for given input prior, prior 
% reinforcement, TFA methods, and target instability. Script relies on
% upper and lower bounds for lambda ranges derived from bStARS to speed 
% computation time (see reference below.)
%% References:
% Miraldi et al. "Leveraging chromatin accessibility data for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
% Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, 
%   R. and Simon, N. -- http://www.stanford.edu/~hastie/glmnet_matlab/
% Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
% Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% INPUTS:
% geneExprMat -- a .mat file containing gene expression data, gene lists
%   (target genes, potential regulators,...), e.g., as generated by
%   importGeneExpGeneLists.m
% tfaMat -- a .mat file containing the prior of TF-gene interactions as
%   well as TFA (prior-based and TF mRNA), e.g., as generated by 
%   integratePrior_estTFA.m
% lambdaBias -- a fractional lambda penalty that a TF-gene interaction in 
%   the prior matrix will have (e.g., penalty term for prior-supported edge
%   is reduced to bias * lambda, where lambda corresponds to the "reference" 
%   penalty applied to edges not in the prior), belongs to range [0,1]
% tfaOpt -- two options recognized:
%   '' --> TFA based on target gene and priors is used, 
%   '_TFmRNA' --> TF based on TF mRNA levels are used
% totSS -- total number of subsamples for the final instability estimates
% targetInstability -- instability cutoff of interest, belongs to range 
%   (0,.5]
% lambdaMin -- initial guess for lambda lower bound containing the lambda
%   that corresponds to the targetInstability, [0, infinity), note glmnet
%   is very slow for lambda < .01
% lambdaMax -- initial guess for lambda lower bound containing the lambda
%   that corresponds to the targetInstability, [0, infinity)
% totLogLambdaSteps -- number of steps per log10 lambda range in
%   bStARS-definied lambda range
% subsampleFrac -- fraction of samples to use for subsampling. Liu et al.,
%   recommend subsample size = floor(100/sqrt(N)), where N = total samples.
%   Given that some TRN inference datasets have < 100 samples, Miraldi et
%   al, used .63*N.
% instabOutMat -- full file name and path for output .mat file, will also
%   be used as file base name for figures and .mat output from bStARS
%   parameter search
% leaveOutSampleList -- a text file, where each line corresponds to a
%   sample condition to be left-out of the inference procedure (e.g., for
%   the purposes of cross-validation)
% bStarsTotSS -- the number of subsamples used to define lambda upper and 
%   lower bound with bStARS.  Muller et al., suggest 2, but a slightly
%   higher number of subsamples (e.g., 5) for bStARS bound derivation can 
%   lead to a tighter lower bound and speed computation overall
% extensionLimit -- in the event that the target instability is outside of
%   user provided [lambdaMin, lambdaMax] range, the range will be extended
%   by an order of magnitude (in the needed direction) up to extensionLimit
%   number of times
%% OUTPUTS:
% instabOutMat -- contains network- and gene-level instabilities,
%   lambdaRange, number of nonzero subsamples per edge (used to rank
%   TF-gene interactions by subsequent scripts)
% instabOutMat.fig + .pdf -- showing network- and gene-level instabilities
%   as a function of final lambda range
% instabOutMat_bStARS -- contains bStARS-derrived upper and lowerbounds for
%   network- and gene-level instabilities
% instabOutMat_bStARS.fig + .pdf -- visualization of upper and lower bounds
%   on target instability as a function of lambda
%% NOTE: Here edges in the prior are treated in a binary manner (present, 
% nonzero, or absent, zero). This code could be modified to include
% real-valued prior edge confidences.

%% load gene expression and TFA
load(geneExprMat)
totSamps = size(targGeneMat,2);
responseMat = targGeneMat;
load(tfaMat)
% have to match prior names with target gene expression and TFA
if tfaOpt
    disp('noTfa option')
    pRegs = pRegsNoTfa;
    pTargs = pTargsNoTfa;
    priorMatrix = priorMatrixNoTfa;
end
[uniNoPriorRegs,uniNoPriorRegInds] = setdiff(potRegs_mRNA,pRegs);
allPredictors = cellstr(strvcat(strvcat(pRegs),strvcat(uniNoPriorRegs)));
totPreds = length(allPredictors);

[vals, targGeneInds, priorGeneInds] = intersect(targGenes,pTargs);
totTargGenes = length(targGenes);
totPRegs = length(strvcat(pRegs));
priorMat = zeros(totTargGenes,totPreds);
priorMat(targGeneInds,1:totPRegs) = priorMatrix(priorGeneInds,:);

%% set input priors and predictors    
predictorMat = [medTfas; potRegMat_mRNA(uniNoPriorRegInds,:)];
if tfaOpt % use the mRNA levels of TFs
    currPredMat = zeros(totPreds,totSamps);
    for prend = 1:totPreds
        prendInd = find(ismember(potRegs_mRNA,allPredictors{prend}));
        currPredMat(prend,:) = potRegMat_mRNA(prendInd,:);
    end        
    predictorMat = currPredMat;
    disp(['TF mRNA used.'])
end
priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*priorMat; %!!! modified from priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*abs(sign(priorMat)); 

if tfaOpt
    %% set lambda penalty to infinity for positive feedback edges where TF 
    % mRNA levels serves both as gene expression and TFA estimate
    for pr = 1:totPreds
        targInd = find(ismember(targGenes,allPredictors{pr})); 
        if length(targInd) % set lambda penalty to infinity, avoid predicting a TF's mRNA based on its own mRNA level
            priorWeightsMat(targInd,pr) = inf; % i.e., target gene is its own predictor
        end
    end    
else % have to set prior inds to zero for TFs in TFA that don't have prior info
    for pr = 1:totPreds        
        if sum(abs(priorMat(:,pr))) == 0 % we have no target edges to estimate TF's TFA
            targInd = find(ismember(targGenes,allPredictors{pr}));
            if length(targInd) % And TF is in the predictor set
                priorWeightsMat(targInd,pr) = inf;
            end
        end
    end
end

is_one = priorWeightsMat == 1;
priorWeightsMat( is_one ) = inf;

%% Check whether to use full gene expression matrix or exclude leave-out set 
if leaveOutSampleList
    disp(['Leave-out set detected: ' leaveOutSampleList])
    % get leave-out set of samples
    fin = fopen(leaveOutSampleList,'r');
    C = textscan(fin,'%s','HeaderLines',0);
    fclose(fin);
    testSamples = C{1};
    testInds = find(ismember(conditionsc,testSamples));
    trainInds = setdiff(1:totSamps,testInds);
else        
    disp(['Full gene expression matrix used.'])
    trainInds = 1:totSamps; % all training samples used
    testInds = [];
end
subsampleSize = floor(subsampleFrac*length(trainInds));

%% bStARS to narrow lambda range around target instability
disp('Estimating lambda bounds for target instability with getMLassoStARSlambdaRangePerGene.m')
figure(1), clf
% bStarsTotSS = 5;  % note bStARS authors recommend 2 subsamples, we reduce search space 
% further by using 5 subsamples at relatively small cost early on
bStarsLogLambdaStep = 10;
lamLog10step = 1/bStarsLogLambdaStep;
logLamRange =  log10(lambdaMin):lamLog10step:log10(lambdaMax);
lambdaRange = [10.^logLamRange];
% get lambda range
[minLambdas, maxLambdas, maxedOut, notSmallEnough,...
    minLambdaNet,maxLambdaNet,maxOutNet,minOutNet,...
    netInstabilitiesLb,netInstabilitiesUb,...
    instabilitiesLb,instabilitiesUb] = ...
    getMLassoStARSlambdaRangePerGene(predictorMat(:,trainInds),responseMat(:,trainInds),...
        priorWeightsMat,lambdaRange,targetInstability,...
        targetInstability,subsampleSize,bStarsTotSS);
% extend lambda range for genes where lambda range was too small
% note extension is limited to "extensionLimit" defined above
needNewRange = maxOutNet + minOutNet;        
extended = 0;
while and(needNewRange,extended < extensionLimit)
    if maxOutNet > 0
        currLambdaMax = lambdaMax*10;
        lambdaMax = currLambdaMax;
    else
        currLambdaMax = lambdaMax;
    end
    if minOutNet > 0
        currLambdaMin = lambdaMin/10;
        lambdaMin = currLambdaMin;
    else
        currLambdaMin = lambdaMin;
    end
    logLamRange =  log10(currLambdaMin):lamLog10step:log10(currLambdaMax);
    lambdaRange = [10.^logLamRange];    
    [minLambdas, maxLambdas, maxedOut, notSmallEnough,...
            minLambdaNet,maxLambdaNet,maxOutNet,minOutNet,...
            netInstabilitiesLb,netInstabilitiesUb,...
            instabilitiesLb,instabilitiesUb] = ...
        getMLassoStARSlambdaRangePerGene(predictorMat(:,trainInds),responseMat(:,trainInds),...
            priorWeightsMat,lambdaRange,targetInstability,...
            targetInstability,subsampleSize,bStarsTotSS);  
    extended = extended + 1;     
    needNewRange = unique([maxedOut; notSmallEnough]);  
    disp(['Extended lambda range ' num2str(extended) ' time(s).'])
end    
figInf = [strrep(instabOutMat,'.mat','') '_bStARS'];
saveas(gcf,figInf,'fig')
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [6 10]);
print('-painters','-dpdf','-r150',[figInf '.pdf'])
disp(figInf)    
save([figInf '.mat'],...
    'minLambdas', 'maxLambdas', 'instabilitiesLb',...
    'instabilitiesUb','netInstabilitiesLb','netInstabilitiesUb',...
    'minLambdaNet','maxLambdaNet','predictorMat','responseMat','priorMat',...
    'lambdaMax','lambdaBias','lambdaRange',...
    'targetInstability','subsampleSize','bStarsTotSS','targGenes','testInds',...
    'priorWeightsMat','allPredictors','conditionsc','trainInds','extended')

disp(['For target instability = ' num2str(targetInstability)])
disp(['Min lambda = ' num2str(minLambdaNet) ', Max lambda = ' num2str(maxLambdaNet)])

%% get instabilities
disp('Estimating instabilities with getMLassoStARSinstabilitiesPerGeneAndNet.m')
lamLog10step = 1/totLogLambdaSteps;
logLamRange =  log10(minLambdaNet):lamLog10step:log10(maxLambdaNet);
lambdaRange = [10.^logLamRange]';

figure(3), clf
% get Network-wise edge instabilities, and gene-wise edge instabilities    
[geneInstabilities,netInstabilities,ssMatrix] = ...
    getMLassoStARSinstabilitiesPerGeneAndNet(predictorMat(:,trainInds),...
    responseMat(:,trainInds),priorWeightsMat,lambdaRange,...
    subsampleSize, totSS);     

figInf = strrep(instabOutMat,'.mat','');
saveas(gcf,figInf,'fig')
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [6 5]);
print('-painters','-dpdf','-r150',[figInf '.pdf'])
disp(figInf)        
save([instabOutMat '.mat'],'-v7.3',...
    'geneInstabilities','netInstabilities','ssMatrix','predictorMat',...
    'responseMat','priorMat','lambdaBias','lambdaRange','trainInds',...
    'subsampleSize','totSS',...
    'targGenes','priorWeightsMat','allPredictors','conditionsc') 
