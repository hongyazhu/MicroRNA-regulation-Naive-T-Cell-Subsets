function vis_topN_TfTfOverlaps(datasetName, tfPairMat, desClusts, topN,...
    titleInf, fullFigOutDir, aveGeneExprMat, annotations, saveFig, ....
    axisFontSize, figureDimensions)
%% vis_topN_TfTfOverlaps(datasetName, tfPairMat, desClusts, topN,...
%     titleInf, fullFigOutDir, aveGeneExprMat, annotations, saveFig, ....
%     axisFontSize, figureDimensions)
%% GOAL: cluster/ visualize the "Top N" most significant TF-TF modules from 
% calc_zscoredTfTfOverlaps.m for a particular clustering solution (find
% total # of clusters using eval_clustSolns_tfTfOverlap.m). 
%% This function produces several plots:
% 3. plot silhouette scores per TF
% 4. similarity matrix normalized overlaps
% 5. similarity matrix based on raw overlaps
% 6. similarity matrix based on raw p-val Fisher exact test
% 7. similarity matrix based on adjusted p-val Fisher exact test
% 8. OPTIONAL heatmap of TF gene expression
% 9 - 8+N. OPTIONAL N heatmaps of enrichments (e.g., p-val matrices of GSEA
% results (TFs X gene sets), as generated by Th17_example/tfTargets_GSEA.sh
% (core TF analysis) or scTRN/tfTargets_GSEA_loop.sh (GSEA of Kegg, GO
% pathways, etc.). Alternatively, GSEA of many sets is likely best visualized using 
% scTRN/visGSEAenrich_heatmaps_comb.m, limiting to the TFs in the output .txt
% BONUS: provides a .txt file of TFs ordered according to dendrogram; first
%   column indicates whether the TF marks the start of a new cluster and
%   the second column corresponds to TF name. This is useful to provide for
%   GSEA visualization code (e.g., scTRN/visGSEAenrich_heatmaps_comb.m)
%% INPUTS:
% datasetName -- a short name for the analysis (e.g., "Th17"), will be used
%   in figure titles and as a base for file names
% tfPairMat -- .mat file of TF-TF normalized overlaps and other statistics, 
%   as calculated by calc_zscoredTfTfOverlaps.m
% desClusts -- desired number of clusters / clustering solution
% topN -- number of "most significant" clusters to visualize
% titleInf -- a string to be included in title figures
% fullFigOutDir -- output directory for figures and .txt file
% aveGeneExprMat -- OPTIONALLY include data to generate a gene expression 
%    heatmap (e.g., as generated by aveGeneExpMatrix_subset_4viz.m), set to
%    empty string '' to omit this figure
% annotations -- OPTIONALLY include a table of gene-set enrichments, as an 
%   N x 3 cell, were each row corresponds to an enrichment analysis:
%   col 1: sign of enriched regulatory edges (1 --> positive, -1 --> inhibition)
%   col 2: nickname for enrichment anlaysis (goes in title & figure file)
%   col 3: tab-delimited table of p-values associated with enrichments (sets X TFs)
%   col 4: cell with the desired ordered of set annotations -- if desired
%       order is unknown, this can be left as an empty string
% saveFig -- binary flag to save figures: 1 --> yes, 0 --> no
% axisFontSize -- heatmap fontsize
% figureDimensions -- 2D vector of (x-size, y-size) dimensions for
%   saving pdfs of the figures, units are inches
%% OUTPUTS:
% Figures corresponding to the 7+ figures described above
% .txt file of clustered TF names and locations, can be used to visualize
%   more enrichments for those TFs, using scTRN/visGSEAenrich_heatmaps_comb.m

%% debugging
% clear all; close all
% restoredefaultpath
% matlabDir = '..';
% cd(fullfile(matlabDir,'Th17_example'))
% addpath(fullfile(matlabDir,'infLassoStARS'))
% addpath(fullfile(matlabDir,'customMatlabFxns'))
% 
% % input inputs:
% inputNetwork = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/ATAC_Th17_bias50_maxComb_cut01.tsv';
% [outDirBase,fileName,ext] = fileparts(inputNetwork);
% tfTargMin = 20;     % only consider TFs with at least this number of targets
% targTfMin = 1;      % only consider targets with at least this number of TFs
% fdrCut = .1;         % cutoff for TF pair inclusion
% edgeOpt = 'comb';   % can be set to one of three options:
% outDir = fullfile(outDirBase,fileName,strjoin({'zOverlaps',...    
%     [edgeOpt 'Edge'],...
%     ['fdr' num2str(100*fdrCut)],...
%     ['tfMin' num2str(tfTargMin)],...
%     ['targMin' num2str(targTfMin)]},'_'));
% xSize = 7;      
% ySize = 12;
% 
% % required inputs:
% datasetName = 'Th17';
% tfPairMat = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/ATAC_Th17_bias50_maxComb_cut01/zOverlaps_combEdge_fdr10_tfMin20_targMin1/tfPair.mat';
% desClusts = 50;     % desired number of clusters 
% topN = 15;
% titleInf = [datasetName ', ' edgeOpt '-edge , FDR: ' num2str(100*fdrCut) '%, min TF, target: ' num2str(tfTargMin) ', ' num2str(targTfMin) ', clust = ' num2str(desClusts)]; % 
% 
% fullFigOutDir = fullfile(outDir,['Top' num2str(topN) '_Figs_clust' num2str(desClusts)]);
% %% OPTIONAL INPUTS:
% aveGeneExprMat = './outputs/processedGeneExpTFA/geneExpHeatmapInputs/Th0_Th17_48hTh.mat';
% % Annotations: each row is a 3-column cell where:
% %   col 1: sign of enriched regulatory edges (1 --> positive, -1 --> inhibition)
% %   col 2: nickname for enrichment anlaysis (goes in title & figure file)
% %   col 3: tab-delimited table of p-values associated with enrichments (sets X TFs)
% %   col 4: cell with the desired ordered of set annotations -- if desired
% %       order is unknown, this can be left as an empty string
% annotations = {... Annotation 1: TFs core due to positive edges
%     1,'posEdgeCore',...
%     'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/GSEA/ATAC_Th17_bias50_maxComb_cut01_Th17set_Praw0p1_dir_wCut0p0_minSet5/Th17set_praw10_up_adjp.txt',...
%     {'Th17 upOnly'; 'Th17 downOnly'};
%     ... Annotation 2: TFs core due to negative edges
%     -1,'negEdgeCore',....
%     'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/GSEA/ATAC_Th17_bias50_maxComb_cut01_Th17set_Praw0p1_dir_wCut0p0_minSet5/Th17set_praw10_down_adjp.txt',...
%     {'Th17 upOnly'; 'Th17 downOnly'};
%     };
% 
% saveFig = 1;        % save figures? 1 --> yes, 0 --> no
% axisFontSize = 4;   % heatmap fontsize
% figureDimensions = [xSize ySize]; % 2D vector of (x-size, y-size) dimensions for
% %   saving pdfs of the figures, units are inches

%% BEGIN function
lineWidth = .33;    % line width for heatmap clusters
padjSat = .001;     % saturation max for heatmap of TF GSEA enrichments
titleFontSize = 14; % title font size
nonSimFactor = 2;   % multiply axisFontSize by this number for non-square 
    % similarity matrix outputs (e.g., dendrogram, silhouette scores,...)

xSize = figureDimensions(1);
ySize = figureDimensions(2);
totEnrich = size(annotations,1);

disp(titleInf)
mkdir(fullFigOutDir)
disp(fullFigOutDir)
figNameBase = fullfile(fullFigOutDir,datasetName);

%% load overlap statistics
load(tfPairMat)
zMatDist = tfPairAnal.zMatDist;
zMat = tfPairAnal.zMat;
sigTfs = tfPairAnal.sigTfs;
totTfs = length(sigTfs);

%% get clustering solns and silhouette scores   
pdis = tril(zMatDist,-1);
pdis = squareform(pdis);
link = linkage(pdis,'ward');    
currSoln = cluster(link,'maxclust',desClusts);
[silScores,aveSilScore,stdSilScore] = evalSilouetteDist(zMatDist,currSoln);

%% cluster within clusters for BEAUTIFUL figures
%% 1. dendrogram
figure(1), clf
subplot(1,4,1:2)
[h, t, horder] = dendrogram(link,0,'Labels',...
    strrep(sigTfs,'_', ' '),'Orientation','left');
set(h,'Color','k')
horderFull = fliplr(horder);
tfPairAnal.horder = horderFull;        
set(gca,'XTick',[],'FontSize',nonSimFactor*axisFontSize)
xlabel(strvcat('   ',['Euclidean Distance']),'Fontsize',nonSimFactor*axisFontSize)
title([titleInf ':: Edges (zOv)'],'FontSize',titleFontSize)
% if saveFig
%     subfigname = [figNameBase '_dend_zOv'];
%     fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); 
%     print('-painters','-dpdf','-r100',[subfigname '.pdf'])
%     saveas(gcf,[subfigname],'fig')
%     disp([subfigname])  
% end

% get start_spots for marking cluster starts and ends in similary matrices    
start_spots = zeros(desClusts,1);
for clust = 1:desClusts
    clustInds = find(currSoln == clust);
    [xx, horderInds] = intersect(horderFull,clustInds);
    start_spots(clust) = min(horderInds);
end    

solnNstats = [];
solnNstats.clustInds = cell(desClusts,1);
solnNstats.clustMemNames = cell(desClusts,1);
solnNstats.clustMedzOv = zeros(desClusts,1);
solnNstats.clustMeanzOv = zeros(desClusts,1);
solnNstats.clustMadPraw = zeros(desClusts,1);
solnNstats.soln = currSoln;
solnNstats.silScores = silScores;
solnNstats.start_spots = start_spots;

zMatDistribution = squareform(zMat);
maxZ = max(zMatDistribution);
totPairs = length(zMatDistribution);
clustPvals = zeros(desClusts,1);
clustSizes = zeros(desClusts,1);
for clust = 1:desClusts
    clustInds = find(currSoln == clust);
    clustSize = length(clustInds);       
    clustSizes(clust) = clustSize;
    solnNstats.clustInds{clust} = clustInds;
    solnNstats.clustMemNames{clust} = strjoin(cellstr(strvcat(upper({sigTfs{clustInds}}))),'_');
%     strjoin(cellstr(strvcat(upper({sigTfs{clustInds}}))),' ')
    currZmat = squareform(zMat(clustInds,clustInds));
    currPairs = length(currZmat);
    currPs = zeros(currPairs,1);

    for pind = 1:currPairs %      currMed = median(currZmat(pind,setdiff(1:clustSize,pind))); % calculate median
        currPs(pind) = length(find(zMatDistribution >= currZmat(pind)))/totPairs;
        length(find(zMatDistribution >= currZmat(pind)))/totPairs;
    end        
    % combine the p-values
    weights = (clustSize/currPairs)*ones(currPairs,1); % these weights with weighted z-method will count each TF vs. TF pairs
    clustPvals(clust) = combineP_stouffersZMethod(currPs,weights); 
%         clustPvals(clust) = 2*combineP_stouffersZMethod(currPs/2,weights); 
%             % heuristic: divide p-vals by 2, then multiple combined p-value
%             % by 2; this gets around the problem of raw empirical p-values
%             % being identically 1 (z-score = infinity and combined p-value
%             % being 1). We get many raw empirical p-values identically 1 
%             % because we set z-scored overlaps that are smaller than
%             % expected by chance to zero, leading to many zero entries in
%             % the data matrix 
end   
solnNstats.clustSizes = clustSizes;
sigTfsOrdered = cellstr(strvcat(upper({sigTfs{horderFull}})));
%     gompers

%% 2. scatter plot: module empirical p-val versus cluster size
% minNzP = min(clustPvals(clustPvals>0));
% minNzPcut = 10.^(floor(log10(minNzP)));
% minNzP = 10.^(floor(log10(minNzP)-1));
% figure(2), clf
[psOrd, inds] = sort(clustPvals,'ascend');
% psOrd(psOrd==0) = minNzP;
% plot(clustSizes(inds),psOrd,'.','Color',[.8 .8 .8],'MarkerSize',50,'LineWidth',3)
% text(clustSizes(inds),psOrd,strrep({solnNstats.clustMemNames{inds}},'_', ' '))
% set(gca,'yscale','log')%,'FontSize',axisFontSize)
% grid on, ax = axis(); grid minor
% axis([ax(1) ax(2)+1 ax(3) ax(4)])
% hold on
% % plot([ax(1) ax(2)+1],minNzPcut*ones(2,1),'Color',.66*[1 1 1],'LineWidth',2)
% xlabel('Cluster Size')%,'FontSize',axisFontSize)
% ylabel('Signficance')%,'FontSize',axisFontSize)
% title(titleInf)%,'FontSize',titleFontSize)
% if saveFig
%     subfigname = [figNameBase '_clustPvals'];
%     fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize xSize]); 
%     print('-painters','-dpdf','-r100',[subfigname '.pdf'])
% %     pause
% %     save2pdf(subfigname,gcf,100)
%     saveas(gcf,[subfigname],'fig')
%     disp([subfigname])  
% end

%% limit subsequent analyses to Top N TF-TF clusters

clusts2keep = inds(1:topN);
topClustInds = [];

clust2remove = setdiff(1:desClusts,clusts2keep);
inds2remove = find(ismember(currSoln(horderFull),clust2remove));
topOrder = horderFull;
topOrder(inds2remove) = [];    
topClustSolns = currSoln(horderFull);
topClustSolns(inds2remove) = [];      
totKeepClusts = length(clusts2keep);
topStarts = zeros(totKeepClusts,1);
for kc = 1:totKeepClusts
    [vals indLoc] = intersect(topOrder,topOrder(find(topClustSolns == clusts2keep(kc))));
    topStarts(kc) = min(indLoc);        
end
totTopTfs = length(topOrder);
topSigTfs = upper(cellstr(strvcat(sigTfs{topOrder})));


%% 3. silhouette scores per TF
figure(3), clf
subplot(1,4,3:4)
barh(flipud(silScores(topOrder)),'FaceColor','k')
axis tight
set(gca,'YTick',1:totTopTfs,'YTickLabel',strvcat(fliplr(topSigTfs)),...
    'FontSize',axisFontSize*nonSimFactor,'TickLength',[0 0])
xlabel('Silhouette Score','FontSize',nonSimFactor*axisFontSize)
title([titleInf],'FontSize',titleFontSize)    
hold on, axis tight, grid on
%     ax = axis();
%     start_spotsS = start_spots;
%     start_spotsS(start_spots == 1) = [];
%     for ss = 1:desClusts        
%         plot(ax(1:2),((totTopTfs-start_spotsS(ss))*[1 1]-.5),'Color',.25*ones(1,3))
%     end    
if saveFig
    subfigname = [figNameBase '_sil'];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); 
    print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,[subfigname],'fig')
    disp([subfigname])
end

%% 4. similarity matrix normalized overlaps
plotMat = tfPairAnal.zMat(topOrder,topOrder)+max(tfPairAnal.zMat(:))*eye(totTopTfs);
currLabels = cellstr(strvcat(topSigTfs));
figure(4), clf
imagesc(plotMat)
if min(plotMat(:)) < 0
    % most recent implementation of calc_zscoredTfTFOverlaps now provides
    % negative values
    colormap redblue
    clim = max(abs(plotMat(:)));
    set(gca,'CLim',[-clim clim]);
else
    colormap hot
    set(gca,'FontSize',axisFontSize,'FontWeight','Bold')
end
axis image
colorbar
hold on
ax = axis();
for ss = 1:topN
    plot(ax(1:2),(topStarts(ss)*[1 1]-.5),'k','LineWidth',lineWidth)
    plot([topStarts(ss)*[1 1]-.5], ax(3:4),'k','LineWidth',lineWidth)
end    
title([titleInf ':: Normalized Overlap'],'FontSize',titleFontSize)        
set(gca,'YTick',1:totTopTfs,'YTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize)%,'FontWeight','Bold')
set(gca,'XTick',1:totTopTfs,'XTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize,'XTickLabelRotation',90)
set(gca,'TickLength',[0 0])
if saveFig
    subfigname = [figNameBase '_zMat'];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,[subfigname],'fig')
    disp([subfigname])  
    % Save ordered TFs
	fid = fopen([subfigname,'_TF_order_top.txt'],'w');
	for ix = 1:length(currLabels)
		if (ismember(ix,topStarts))
			currStr = [num2str(1),'	',char(currLabels(ix))];
		else
			currStr = [num2str(0),'	',char(currLabels(ix))];
		end
		fprintf(fid,'%s\n',currStr);
	end
	fclose(fid);
end

%% 5. similarity matrix based on raw overlaps
plotMat = tfPairAnal.sigPerOverlaps(topOrder,topOrder)+max(tfPairAnal.sigPerOverlaps(:))*eye(totTopTfs);
currLabels = cellstr(topSigTfs);
figure(5), clf
imagesc(plotMat)
colormap hot
set(gca,'FontSize',axisFontSize,'FontWeight','Bold')
axis image
colorbar
title([titleInf ':: Raw Overlaps'],'FontSize',titleFontSize)        
set(gca,'YTick',1:totTopTfs,'YTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize)%,'FontWeight','Bold')
set(gca,'XTick',1:totTopTfs,'XTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize,'XTickLabelRotation',90)
set(gca,'TickLength',[0 0])    
hold on
ax = axis();
for ss = 1:topN
    plot(ax(1:2),(topStarts(ss)*[1 1]-.5),'r','LineWidth',lineWidth)
    plot([topStarts(ss)*[1 1]-.5], ax(3:4),'r','LineWidth',lineWidth)
end    

if saveFig
    subfigname = [figNameBase '_ov'];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,[subfigname],'fig')
    disp([subfigname])  
end


%% 6. similarity matrix based on raw p-val Fisher exact test
plotMat = (-tfPairAnal.sigRawSigs(topOrder,topOrder))+max(-tfPairAnal.sigRawSigs(:))*eye(totTopTfs);
currLabels = topSigTfs;
figure(6), clf
imagesc(plotMat)
colormap(flipud(gray))
set(gca,'FontSize',axisFontSize,'FontWeight','Bold')
axis image
colorbar
hold on
ax = axis();
for ss = 1:topN
    plot(ax(1:2),(topStarts(ss)*[1 1]-.5),'k','LineWidth',lineWidth)
    plot([topStarts(ss)*[1 1]-.5], ax(3:4),'k','LineWidth',lineWidth)
end
title([titleInf ':: -log_{10}(P_{raw})'],'FontSize',titleFontSize)
set(gca,'YTick',1:totTopTfs,'YTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize)%,'FontWeight','Bold')
set(gca,'XTick',1:totTopTfs,'XTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize,'XTickLabelRotation',90)
set(gca,'TickLength',[0 0])

if saveFig
    subfigname = [figNameBase '_Praw'];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,[subfigname],'fig')
    disp([subfigname])  
end

%% 7. similarity matrix based on adjusted p-val Fisher exact test
plotMat = (tfPairAnal.sigSigs(topOrder,topOrder)+max(tfPairAnal.sigSigs(:))*eye(totTopTfs));
currLabels = topSigTfs;
figure(7), clf
imagesc(plotMat)
colormap(flipud(gray))
set(gca,'FontSize',axisFontSize,'FontWeight','Bold')
axis image
colorbar

hold on
ax = axis();
for ss = 1:topN
    plot(ax(1:2),(topStarts(ss)*[1 1]-.5),'k','LineWidth',lineWidth)
    plot([topStarts(ss)*[1 1]-.5], ax(3:4),'k','LineWidth',lineWidth)
end    
title([titleInf ':: (-log_{10}(P_{adj}))'],'FontSize',titleFontSize)        
set(gca,'YTick',1:totTopTfs,'YTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize)%,'FontWeight','Bold')
set(gca,'XTick',1:totTopTfs,'XTickLabel',strvcat(currLabels),...
    'FontSize',axisFontSize,'XTickLabelRotation',90)
set(gca,'TickLength',[0 0])

if saveFig
    subfigname = [figNameBase '_pAdj'];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,[subfigname],'fig')
    disp([subfigname])  
end    

%% 8. OPTIONAL heatmap of TF gene expression
if length(aveGeneExprMat) > 0
    load(aveGeneExprMat)
    [totGenes, totConds] = size(zAveCounts);        
    geneExprOrd = zeros(totTopTfs,totConds);
    [names,porder,corder] = intersect(upper(genesc),topSigTfs);
    transVals(corder,:) = zAveCounts(porder,:);

    figure(8), clf
    subplot(1,3,2:3)
    imagesc(transVals')
    colormap(redblue)
    set(gca,'XTick',1:totTopTfs,'XTickLabel',topSigTfs,...
        'FontSize',axisFontSize,'FontWeight','Bold','XTickLabelRotation',90)
    set(gca,'YTick',1:totConds,'YTickLabel',strrep(condPrintNames,'_',' '),...
        'FontSize',axisFontSize)
    axis image
    clim = max(abs(transVals(:)));
    set(gca,'CLim',[-clim clim],'TickLength',[0 0])
    colorbar('YTick',[-clim 0 clim],'YTickLabel',{['<-' roundstring1(clim)];'0';['>' roundstring1(clim)]},...
        'FontWeight','Bold')
    title('Z-scored Gene Expression','FontSize',titleFontSize)
    hold on
    ax = axis();
    for ss = 1:topN
        plot([topStarts(ss)*[1 1]-.5], ax(3:4),'k','LineWidth',lineWidth)
    end
    for ss = find(startSpotsConds)
        plot(ax(1:2),(ss*[1 1]-.5),'k','LineWidth',lineWidth)
    end    

    if saveFig
        subfigname = [figNameBase '_expZ'];
        fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [ySize xSize]);  % x / y reversed since this is flipped relative to other heatmaps
        print('-painters','-dpdf','-r100',[subfigname '.pdf'])
        saveas(gcf,[subfigname],'fig')
        disp([subfigname])  
    end 
end

%% make heatmaps of gene set enrichments per TF
%% 9 - 8+N. OPTIONAL N heatmaps of enrichments (e.g., p-val matrices of GSEA
% results (TFs X gene sets), as generated by Th17_example/tfTargets_GSEA.sh
% (core TF analysis) or scTRN/tfTargets_GSEA_loop.sh (GSEA of Kegg, GO
% pathways, etc.); NOTE: postive and negative enrichments assumed
for gs = 1:totEnrich
    %% load file and get in clustering-solution order
    fcDir = annotations{gs,1};
    setName = annotations{gs,2};
    pvalFile = annotations{gs,3};
    orderedSets = annotations{gs,4};
    
    fid = fopen(pvalFile);  % get first line, TFs
    tline=fgetl(fid);
    tfTmps = strvcat(strsplit(tline,'\t'));
    fclose(fid);
    tfsTmpc = upper(cellstr(tfTmps));
    totTfsTmp = length(tfsTmpc);
    fid = fopen(pvalFile);  % get pvals + sets
    C = textscan(fid,['%s' repmat('%f',1,totTfsTmp)],'Delimiter','\t','Headerlines',1);
    setsTmp = strrep((C{1}),'_', ' ');
    
    padjsTmp = [C{2:end}];
    
    % optionally reorder and filter sets
    if isempty(orderedSets)
        padjs = padjTmp;
        sets = setsTmp;
    else
        totIncSets = size(orderedSets,1);
        padjs = zeros(totIncSets,totTfsTmp);
        for sind = 1:totIncSets
            keepInd = find(ismember(setsTmp,orderedSets{sind}));
            padjs(sind,:) = padjsTmp(keepInd,:);
        end
        sets = orderedSets;
    end
    totSets = length(sets);
    
    % reorder TFs
    padjsOrd = ones(totTopTfs,totSets);
    [names,porder,corder] = intersect(tfsTmpc,topSigTfs);
    padjsOrd(corder,:) = padjs(:,porder)';
    transVals = fcDir * -log10(max(padjsOrd,padjSat));        
    
    figure(8+gs), clf
    imagesc(transVals)
    colormap redblue
    set(gca,'YTick',1:totTopTfs,'YTickLabel',topSigTfs,...
        'FontSize',nonSimFactor*axisFontSize,'FontWeight','Bold')
    set(gca,'XTick',1:totSets,'XTickLabel',strrep( sets,'_', ' '),...
        'XTickLabelRotation',90,'FontSize',nonSimFactor*axisFontSize)
    axis image
    set(gca,'CLim',[log10(padjSat) -log10(padjSat)],'TickLength',[0 0])
    try
        colorbar('YTick',fcDir*[0 1 -log10(padjSat)],'YTickLabel',{'P_{adj} = 1';'P_{adj} = .1';['P_{adj} <= ' roundstring2(padjSat)]},...
            'FontWeight','Bold')
    catch
        colorbar('YTick',fliplr(fcDir*[0 1 -log10(padjSat)]),'YTickLabel',flipud({'P_{adj} = 1';'P_{adj} = .1';['P_{adj} <= ' roundstring2(padjSat)]}),...
            'FontWeight','Bold')
    end
    title(setName,'FontSize',titleFontSize)
    hold on
    ax = axis();
    for ss = 1:topN
        plot(ax(1:2),(topStarts(ss)*[1 1]-.5),'k','LineWidth',lineWidth)
    end    

    if saveFig
        subfigname = [figNameBase '_' setName];
        fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); 
        print('-painters','-dpdf','-r100',[subfigname '.pdf'])
        saveas(gcf,[subfigname],'fig')
        disp([subfigname])  
    end     

end

subfigname = [figNameBase '_orderedTFs.txt'];
fout = fopen(subfigname,'w');
for tind = 1:length(sigTfsOrdered)
    if intersect(start_spots,tind)> 0
        fprintf(fout,['1\t' sigTfsOrdered{tind} '\n']);
    else
        fprintf(fout,['0\t' sigTfsOrdered{tind} '\n']);
    end
end
fclose(fout);
disp(subfigname)
