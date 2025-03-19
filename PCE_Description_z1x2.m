function clusterDescr = ...
    PCE_Description_z1x2(DimStruct, clusterMatrix_neg, clusterMatrix_pos, clusterMeasure_neg, clusterMeasure_pos, clustSignThreshold, clustMaxMeasure_H0, statMatrix, eSizMatrix)

% PCE_Description_z1x2z1 describes the clusters 
%
% INPUT: 
%
% DimStruct         z1=Chan, x2=Time
%                   For example...
%                     DimStruct.z1_lbl   = 'Channel';
%                     DimStruct.z1_contFlag = 0;
%                     DimStruct.z1_chanlocs =; chanlocs structure in the eeglab format
%                     DimStruct.x2_lbl   = 'Time';
%                     DimStruct.x2_units = 's';
%                     DimStruct.x2_vec   = activity_time;

%
% OUTPUT:
%
% 
% Author: Germano Gallicchio 
% germano.gallicchio@gmail.com

%% debugging cell

% clustMaxMass_H0 = clustMaxMass1_H0;


%% sanity checks

% size of pos and neg clusterMatrix must be the same 
if ~isequal(size(clusterMatrix_neg), size(clusterMatrix_pos))
    error('size of pos and neg clusterMatrix must be the same')
end

% size of statMatrix and eSizMatrix must be the same 
if ~isequal(size(statMatrix), size(eSizMatrix))
    error('size of stat and eSiz must be the same')
end

%% get data

if DimStruct.z1_contFlag==0
    %disp('z1 is not continuous (e.g., EEG channels)')
    chanlocs = DimStruct.z1_chanlocs;
else
    error('z1 is continuous (e.g., time, freq). If you have two continuous variables, use the y1x2 design')
end

nz1 = length(DimStruct.z1_chanlocs);
nx2 = length(DimStruct.x2_vec);


[z1Matrix,     x2Matrix]       = ndgrid(1:nz1,1:nx2);
[z1Matrix_lbl, x2Matrix_units] = ndgrid(string({chanlocs.labels}),DimStruct.x2_vec);

%%

% cluster measure (e.g., size, mass, rob mass)
clusterMeasure = [clusterMeasure_neg clusterMeasure_pos];
nClust = length(clusterMeasure);

% cluster ID
% each cluster is identified by an ordinal number (negative for negative clusters)
% the number itself is not meaningful. but it serves as identifier.
% negative
clusterID_neg = nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    clusterID_neg(clIdx) = clIdx*-1;
end
% positive
clusterID_pos = nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    clusterID_pos(clIdx) = clIdx;
end
clusterID = [clusterID_neg clusterID_pos];

% pvalues at cluster level (this corresponds with the inferential test)
% p value for each cluster (how many more extreme values are there in the H0 dist?)
% negative
pval_neg = nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    % % BUG. DELETE THIS LINE pval_neg(clIdx) = sum(clustMaxMeasure_H0>abs(clusterMeasure_neg(clIdx))) / length(clustMaxMeasure_H0);
    pval_neg(clIdx) = sum(abs(clusterMeasure_neg(clIdx))<clustMaxMeasure_H0) / length(clustMaxMeasure_H0);
end
% positive
pval_pos = nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    pval_pos(clIdx) = sum(clusterMeasure_pos(clIdx)<clustMaxMeasure_H0) / length(clustMaxMeasure_H0);
end
% concatenate negative and positive
pval = [pval_neg pval_pos];


% x2 information
% find 0, 25, 50, 75, 100th prctile of the x2 dimension (e.g., time) for each cluster
% negative
x2_0_neg  = nan(1,length(clusterMeasure_neg));
x2_25_neg = nan(1,length(clusterMeasure_neg));
x2_50_neg = nan(1,length(clusterMeasure_neg));
x2_75_neg = nan(1,length(clusterMeasure_neg));
x2_100_neg= nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    x2_0_neg(clIdx)  = prctile(x2Matrix_units(idx),0);
    x2_25_neg(clIdx) = prctile(x2Matrix_units(idx),25);
    x2_50_neg(clIdx) = median(x2Matrix_units(idx));
    x2_75_neg(clIdx) = prctile(x2Matrix_units(idx),75);
    x2_100_neg(clIdx)= prctile(x2Matrix_units(idx),100);
end
% positive
x2_0_pos  = nan(1,length(clusterMeasure_pos));
x2_25_pos = nan(1,length(clusterMeasure_pos));
x2_50_pos = nan(1,length(clusterMeasure_pos));
x2_75_pos = nan(1,length(clusterMeasure_pos));
x2_100_pos= nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    x2_0_pos(clIdx)  = prctile(x2Matrix_units(idx),0);
    x2_25_pos(clIdx) = prctile(x2Matrix_units(idx),25);
    x2_50_pos(clIdx) = median(x2Matrix_units(idx));
    x2_75_pos(clIdx) = prctile(x2Matrix_units(idx),75);
    x2_100_pos(clIdx)= prctile(x2Matrix_units(idx),100);
end
% concatenate negative and positive
x2_0  = [ x2_0_neg    x2_0_pos ];
x2_25 = [ x2_25_neg   x2_25_pos ];
x2_50 = [ x2_50_neg   x2_50_pos ];
x2_75 = [ x2_75_neg   x2_75_pos ];
x2_100= [ x2_100_neg  x2_100_pos ];



% z1, channel information
% for each cluster, list channels representing all occurrences
% negative
chanInCluster_neg = cell(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    chanInCluster_neg{clIdx} = sort(unique(z1Matrix_lbl(idx)))';
end
% positive
chanInCluster_pos = cell(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    chanInCluster_pos{clIdx} = sort(unique(z1Matrix_lbl(idx)))';
end
% concatenate negative and positive
chanInCluster = [chanInCluster_neg chanInCluster_pos];
chanInClusterNum = cellfun(@length,chanInCluster);



% effect size (e.g., Cohen's d, correlation coefficient)
% eSiz within each cluster (this does not correspond with any inferential test, purely descriptive)
%   - Mdn (median) 
%   - Ext (most extreme) of eSiz across all pixels within each cluster
%   - weighted Medoid
% negative
eSizMdn_neg = nan(1,length(clusterMeasure_neg));
eSizExt_neg = nan(1,length(clusterMeasure_neg));
eSizExt_z1_lbl_neg   = repmat("",1,length(clusterMeasure_neg));
eSizExt_z1_idx_neg   = nan(1,length(clusterMeasure_neg));
eSizExt_x2_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizExt_x2_units_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_z1_lbl_neg   = repmat("",1,length(clusterMeasure_neg));
eSizMedoid_z1_idx_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_x2_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_x2_units_neg = nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    % median
    eSizMdn_neg(clIdx) = median(eSizMatrix(idx));
    % most extreme
    [eSizExt_neg(clIdx),   ExtIdx] = min(eSizMatrix(idx));
    temp = x2Matrix(idx);           eSizExt_x2_pnt_neg(clIdx)   = temp(ExtIdx);
    temp = x2Matrix_units(idx);     eSizExt_x2_units_neg(clIdx) = temp(ExtIdx);
    temp = z1Matrix_lbl(idx);       eSizExt_z1_lbl_neg(clIdx)   = temp(ExtIdx);
    temp = z1Matrix(idx);           eSizExt_z1_idx_neg(clIdx)   = temp(ExtIdx);
    % weighted medoid 
    % step 1. weighted medoid to find representative channel (the channel that minimizes the dissimilarity between channels, with dissimilarity = binned angular distance, and with bias towards channels with largest cluster mass portion)
    statMat = eSizMatrix;
    clustIdxMat = idx;
    angularDistMat = DimStruct.z1_D;
    [~, z1_wMedIdx] = PCE_z1Medoids_z1x2(statMat, clustIdxMat,angularDistMat);
    % step 2. consider only points within the representative channel, then find weighted medoid as representative point (the point that minimizes the dissimilarity, with dissimilarity Euclidean distance in the remaining dimensions
    % subsetting idx and eSizMatrix so that only the representative channel is considered
    idx_slice        = squeeze(idx(z1_wMedIdx,:));   
    eSizMatrix_slice = squeeze(eSizMatrix(z1_wMedIdx,:));   
    [~, x2_wMedIdx] = PCE_x2Medoids(eSizMatrix_slice,idx_slice);
    eSizMedoid_z1_idx_neg(clIdx)   = z1_wMedIdx;
    eSizMedoid_x2_pnt_neg(clIdx)   = x2_wMedIdx;
    eSizMedoid_z1_lbl_neg(clIdx)   = z1Matrix_lbl(z1_wMedIdx,x2_wMedIdx);
    eSizMedoid_x2_units_neg(clIdx) = x2Matrix_units(z1_wMedIdx,x2_wMedIdx);
    eSizMedoid_neg(clIdx)  = eSizMatrix(z1_wMedIdx,x2_wMedIdx);
end
% positive
eSizMdn_pos = nan(1,length(clusterMeasure_pos));
eSizExt_pos = nan(1,length(clusterMeasure_pos));
eSizExt_z1_lbl_pos   = repmat("",1,length(clusterMeasure_pos));
eSizExt_z1_idx_pos   = nan(1,length(clusterMeasure_pos));
eSizExt_x2_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizExt_x2_units_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_z1_lbl_pos   = repmat("",1,length(clusterMeasure_pos));
eSizMedoid_z1_idx_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_x2_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_x2_units_pos = nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    % median
    eSizMdn_pos(clIdx) = median(eSizMatrix(idx));
    % most extreme
    [eSizExt_pos(clIdx),   ExtIdx] = max(eSizMatrix(idx));
    temp = z1Matrix_lbl(idx);       eSizExt_z1_lbl_pos(clIdx)   = temp(ExtIdx);
    temp = z1Matrix(idx);           eSizExt_z1_idx_pos(clIdx)   = temp(ExtIdx);
    temp = x2Matrix(idx);           eSizExt_x2_pnt_pos(clIdx)   = temp(ExtIdx);
    temp = x2Matrix_units(idx);     eSizExt_x2_units_pos(clIdx) = temp(ExtIdx);
    % weighted medoid 
    % step 1. weighted medoid to find representative channel (the channel that minimizes the dissimilarity between channels, with dissimilarity = binned angular distance, and with bias towards channels with largest cluster mass portion)
    statMat = eSizMatrix;
    clustIdxMat = idx;
    angularDistMat = DimStruct.z1_D;
    [~, z1_wMedIdx] = PCE_z1Medoids_z1x2(statMat, clustIdxMat,angularDistMat);
    % step 2. consider only points within the representative channel, then find weighted medoid as representative point (the point that minimizes the dissimilarity, with dissimilarity Euclidean distance in the remaining dimensions
    % subsetting idx and eSizMatrix so that only the representative channel is considered
    idx_slice        = squeeze(idx(z1_wMedIdx,:));   
    eSizMatrix_slice = squeeze(eSizMatrix(z1_wMedIdx,:));   
    [~, x2_wMedIdx] = PCE_x2Medoids(eSizMatrix_slice,idx_slice);
    eSizMedoid_z1_idx_pos(clIdx)   = z1_wMedIdx;
    eSizMedoid_x2_pnt_pos(clIdx)   = x2_wMedIdx;
    eSizMedoid_z1_lbl_pos(clIdx)   = z1Matrix_lbl(z1_wMedIdx,x2_wMedIdx);
    eSizMedoid_x2_units_pos(clIdx) = x2Matrix_units(z1_wMedIdx,x2_wMedIdx);
    eSizMedoid_pos(clIdx)  = eSizMatrix(z1_wMedIdx,x2_wMedIdx);
end
% concatenate negative and positive
eSizMdn = [ eSizMdn_neg   eSizMdn_pos ];
eSizExt = [ eSizExt_neg   eSizExt_pos ];
eSizExt_z1_lbl   = [ eSizExt_z1_lbl_neg    eSizExt_z1_lbl_pos   ];
eSizExt_z1_idx   = [ eSizExt_z1_idx_neg    eSizExt_z1_idx_pos   ];
eSizExt_x2_pnt   = [ eSizExt_x2_pnt_neg    eSizExt_x2_pnt_pos   ];
eSizExt_x2_units = [ eSizExt_x2_units_neg  eSizExt_x2_units_pos ];
eSizMedoid = [ eSizMedoid_neg   eSizMedoid_pos ];
eSizMedoid_z1_lbl   = [ eSizMedoid_z1_lbl_neg    eSizMedoid_z1_lbl_pos  ];
eSizMedoid_z1_idx   = [ eSizMedoid_z1_idx_neg    eSizMedoid_z1_idx_pos  ];
eSizMedoid_x2_pnt   = [ eSizMedoid_x2_pnt_neg    eSizMedoid_x2_pnt_pos  ];
eSizMedoid_x2_units = [ eSizMedoid_x2_units_neg  eSizMedoid_x2_units_pos];





% stat values at cluster level (this does not correspond with any test, purely descriptive)
% Mdn (median) and Ext (most extreme) of stat (e.g., tvalue, corr coefficients) across all pixels within each cluster
% negative
statMdn_neg = nan(1,length(clusterMeasure_neg));
statExt_neg = nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    statMdn_neg(clIdx) = median(statMatrix(idx));
    statExt_neg(clIdx) = min(statMatrix(idx));
end
% positive
statMdn_pos = nan(1,length(clusterMeasure_pos));
statExt_pos = nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    statMdn_pos(clIdx) = median(statMatrix(idx));
    statExt_pos(clIdx) = max(statMatrix(idx));
end
% concatenate negative and positive
statMdn = [ statMdn_neg   statMdn_pos ];
statExt = [ statExt_neg   statExt_pos ];




clusterDescr = struct();
for clIdx = 1:nClust
    clusterDescr(clIdx).id        = clusterID(clIdx);
    clusterDescr(clIdx).measure   = clusterMeasure(clIdx);
    clusterDescr(clIdx).threshold = clustSignThreshold;
    clusterDescr(clIdx).pvalue = pval(clIdx);
    
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_0'  ]) = x2_0(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_25' ]) = x2_25(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_50' ]) = x2_50(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_75' ]) = x2_75(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_100']) = x2_100(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_range']) = ['[' num2str(x2_0(clIdx),'%.2f') ', ' num2str(x2_100(clIdx),'%.2f') ']'];

    clusterDescr(clIdx).([DimStruct.z1_lbl '_count' ]) = chanInClusterNum(clIdx);
    clusterDescr(clIdx).([DimStruct.z1_lbl '_lbl' ])   = strjoin(chanInCluster{clIdx});

    clusterDescr(clIdx).eSiz_Mdn = eSizMdn(clIdx);

    clusterDescr(clIdx).eSiz_Ext = eSizExt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizExt_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizExt_x2_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.z1_lbl '_' 'lbl'])               = eSizExt_z1_lbl(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.z1_lbl '_' 'idx'])               = eSizExt_z1_idx(clIdx);

    clusterDescr(clIdx).eSiz_wMedoid = eSizMedoid(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizMedoid_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizMedoid_x2_units(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.z1_lbl '_' 'lbl'])               = eSizMedoid_z1_lbl(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.z1_lbl '_' 'idx'])               = eSizMedoid_z1_idx(clIdx);

    clusterDescr(clIdx).stat_Mdn = statMdn(clIdx);
    clusterDescr(clIdx).stat_Ext = statExt(clIdx);
end

clusterDescr = struct2table(clusterDescr);
if nClust>0
    [~, sortIdx] = sort(abs(clusterDescr.measure),'descend');
    clusterDescr = clusterDescr(sortIdx,:);
end
