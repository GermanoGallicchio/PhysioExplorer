function clusterDescr = ...
    PE_Description_y1x2z3( ...
    DimStruct, ...
    PE_parameters, ...
    clustIDList, ...
    clusterMatrix, ...
    clusterMeasure_obs, ...
    clustThreshold, ...
    clusterMetricsH0, ...
    statMatrix)

% PE_Description_y1x2z3 describes the clusters 
%
% INPUT: 
%
%   DimStruct   Structure describing the dimensions of the data,
%               with fields as in previous scripts
%
%   PE_parameters        Structure containing analysis parameters, including:
%       .clusterMetricChoice (string) Name of cluster metric to use for sorting
%
%   clustIDList          Vector of cluster IDs to describe (numeric)
%   clusterMatrix        3D matrix labeling each (y1,x2,z3) point with cluster ID (same size as data)
%   clusterMeasure_obs   Struct with observed values for each cluster-level measure (fields: e.g., 'size', 'mass')
%   clustThreshold       Struct with threshold values for each cluster-level measure (fields: e.g., 'size_oneTail')
%   clusterMetricsH0     Struct with null distributions for each cluster-level measure (fields: arrays of H0 values)
%   statMatrix           3D matrix of statistical values (e.g., t-values, correlation coefficients) (same size as data)
%
% OUTPUT:
%
%   clusterDescr         Table summarizing descriptive statistics for each cluster.
%                       Each row represents one cluster and contains:
%                        - Cluster ID and cluster-level measure values
%                        - p-values and thresholds for each measure
%                        - Range and quantiles for continuous dimensions y1 and x2
%                        - Number and names of involved channels in dimension z3
%                        - Median and extreme values of effect size/statistics within the cluster
%                        - Localization info for extreme points and weighted medoids
%                        - All fields are labeled using the provided DimStruct labels/units
%
%
% 
% Author:   Germano Gallicchio (germano.gallicchio@gmail.com)
%

%% input validaiton / sanity checks

% size of statMatrix and eSizMatrix must be the same 
% if ~isequal(size(statMatrix), size(eSizMatrix))
%     error('size of stat and eSiz must be the same')
% end

%%
% TO DO: add to PE_description to only use statMatrix, whatever it is
% TO DO: tidy up this function a bit
%   - prctiles for y1 and x2 can be computed in the same loop

%% get data

ny1 = length(DimStruct.y1_vec);
nx2 = length(DimStruct.x2_vec);
nz3 = length(DimStruct.z3_chanlocs);

[y1Matrix,       x2Matrix,       z3Matrix,    ] = ndgrid(1:ny1,1:nx2,1:nz3);
[y1Matrix_units, x2Matrix_units, z3Matrix_lbl ] = ndgrid(DimStruct.y1_vec,DimStruct.x2_vec,string({DimStruct.z3_chanlocs.labels}));

%%

% cluster measure (e.g., size, mass)
nClust = length(clustIDList);

% oneTail pvalues at cluster level (this corresponds with the inferential test)
% p value for each cluster (how many more extreme values are there in the H0 dist?)
clusterMeasure_lbl = fieldnames(clusterMeasure_obs);
clusterMeasure_lbl = clusterMeasure_lbl(~strcmp(clusterMeasure_lbl,'id')); % remove the 'id' column
clusterMeasure_num = length(clusterMeasure_lbl);
pval = nan(nClust,clusterMeasure_num);
for msIdx = 1:clusterMeasure_num
    H0distribution = [clusterMetricsH0.(clusterMeasure_lbl{msIdx})]; 
    obsVal = clusterMeasure_obs.(clusterMeasure_lbl{msIdx}); % observed cluster's measure
    for clIdx = 1:nClust
        pval(clIdx,msIdx) = sum(abs(H0distribution)>=abs(obsVal(clIdx))) / length(H0distribution);  
    end
end
%%

% y1 information
% find 0, 25, 50, 75, 100th prctile of the y1 dimension (e.g., freq) for each cluster
% negative

y1_0  = nan(1,nClust);
y1_25 = nan(1,nClust);
y1_50 = nan(1,nClust);
y1_75 = nan(1,nClust);
y1_100= nan(1,nClust);
for clIdx = 1:nClust
    idx = clusterMatrix==clustIDList(clIdx);
    y1_0(clIdx)  = prctile(y1Matrix_units(idx),0);
    y1_25(clIdx) = prctile(y1Matrix_units(idx),25);
    y1_50(clIdx) = median(y1Matrix_units(idx));
    y1_75(clIdx) = prctile(y1Matrix_units(idx),75);
    y1_100(clIdx)= prctile(y1Matrix_units(idx),100);
end

% x2 information
% find 0, 25, 50, 75, 100th prctile of the x2 dimension (e.g., time) for each cluster
% negative
x2_0   = nan(1,nClust);
x2_25  = nan(1,nClust);
x2_50  = nan(1,nClust);
x2_75  = nan(1,nClust);
x2_100 = nan(1,nClust);
for clIdx = 1:nClust
    idx = clusterMatrix==clustIDList(clIdx);
    x2_0(clIdx)  = prctile(x2Matrix_units(idx),0);
    x2_25(clIdx) = prctile(x2Matrix_units(idx),25);
    x2_50(clIdx) = median(x2Matrix_units(idx));
    x2_75(clIdx) = prctile(x2Matrix_units(idx),75);
    x2_100(clIdx)= prctile(x2Matrix_units(idx),100);
end



% z3, channel information
% for each cluster, list channels representing all occurrences
chanInCluster = cell(1,nClust);
for clIdx = 1:nClust
    idx = clusterMatrix==clustIDList(clIdx);
    chanInCluster{clIdx} = unique(z3Matrix_lbl(idx))';
end
chanInClusterNum = cellfun(@length,chanInCluster);



% effect size (e.g., Cohen's d, correlation coefficient)
% eSiz within each cluster (this does not correspond with any inferential test, purely descriptive)
%   - Mdn (median) 
%   - Ext (most extreme) of eSiz across all pixels within each cluster
%   - weighted Medoid
eSizMdn = nan(1,nClust);
eSizExt = nan(1,nClust);
eSizExt_y1_pnt   = nan(1,nClust);
eSizExt_y1_units = nan(1,nClust);
eSizExt_x2_pnt   = nan(1,nClust);
eSizExt_x2_units = nan(1,nClust);
eSizExt_z3_lbl   = repmat("",1,nClust);
eSizExt_z3_idx   = nan(1,nClust);
eSizMedoid = nan(1,nClust);
eSizMedoid_y1_pnt   = nan(1,nClust);
eSizMedoid_y1_units = nan(1,nClust);
eSizMedoid_x2_pnt   = nan(1,nClust);
eSizMedoid_x2_units = nan(1,nClust);
eSizMedoid_z3_lbl   = repmat("",1,nClust);
eSizMedoid_z3_idx   = nan(1,nClust);
for clIdx = 1:nClust
    idx = clusterMatrix==clustIDList(clIdx);
    % median
    eSizMdn(clIdx) = median(statMatrix(idx));
    % most extreme
    if sign(clustIDList(clIdx)) < 0
        [eSizExt(clIdx),   ExtIdx] = min(statMatrix(idx)); % min for negative clusters
    else
        [eSizExt(clIdx),   ExtIdx] = max(statMatrix(idx)); % max for positive clusters
    end
    temp = y1Matrix(idx);           eSizExt_y1_pnt(clIdx)   = temp(ExtIdx);
    temp = y1Matrix_units(idx);     eSizExt_y1_units(clIdx) = temp(ExtIdx);
    temp = x2Matrix(idx);           eSizExt_x2_pnt(clIdx)   = temp(ExtIdx);
    temp = x2Matrix_units(idx);     eSizExt_x2_units(clIdx) = temp(ExtIdx);
    temp = z3Matrix_lbl(idx);       eSizExt_z3_lbl(clIdx)   = temp(ExtIdx);
    temp = z3Matrix(idx);           eSizExt_z3_idx(clIdx)   = temp(ExtIdx);
    % weighted medoid 
    % step 1. weighted medoid to find representative channel (the channel that minimizes the dissimilarity between channels, with dissimilarity = binned angular distance, and with bias towards channels with largest cluster mass portion)
    statMat = statMatrix;
    clustIdxMat = idx;
    angularDistMat = DimStruct.z3_D;
    [~, z3_wMedIdx] = PE_z3Medoids_y1x2z3(statMat, clustIdxMat,angularDistMat);
    % step 2. consider only points within the representative channel, then find weighted medoid as representative point (the point that minimizes the dissimilarity, with dissimilarity Euclidean distance in the remaining dimensions
    % subsetting idx and eSizMatrix so that only the representative channel is considered
    statMat     = statMatrix(:,:,z3_wMedIdx);   
    clustIdxMat = idx(:,:,z3_wMedIdx);   
    [~, ~, y1_wMedIdx, x2_wMedIdx] = PE_y1x2Medoids(statMat,clustIdxMat);
    eSizMedoid_y1_pnt(clIdx)   = y1_wMedIdx;
    eSizMedoid_x2_pnt(clIdx)   = x2_wMedIdx;
    eSizMedoid_z3_idx(clIdx)   = z3_wMedIdx;
    eSizMedoid_y1_units(clIdx) = y1Matrix_units(y1_wMedIdx,x2_wMedIdx,z3_wMedIdx);
    eSizMedoid_x2_units(clIdx) = x2Matrix_units(y1_wMedIdx,x2_wMedIdx,z3_wMedIdx);
    eSizMedoid_z3_lbl(clIdx)   = z3Matrix_lbl(y1_wMedIdx,x2_wMedIdx,z3_wMedIdx);
    eSizMedoid(clIdx)  = statMatrix(y1_wMedIdx,x2_wMedIdx,z3_wMedIdx);
end


% stat values at cluster level (this does not correspond with any test, purely descriptive)
% Mdn (median) and Ext (most extreme) of stat (e.g., tvalue, corr coefficients) across all pixels within each cluster
statMdn = nan(1,nClust);
statExt = nan(1,nClust);
for clIdx = 1:nClust
    idx = clusterMatrix==clustIDList(clIdx);
    statMdn(clIdx) = median(statMatrix(idx));
    if sign(clustIDList(clIdx)) < 0
        statExt(clIdx) = min(statMatrix(idx));  % min for negative
    else
        statExt(clIdx) = max(statMatrix(idx));  % max for positive
    end
    
    
end



clusterDescr = struct();
for clIdx = 1:nClust
    clusterDescr(clIdx).id        = clustIDList(clIdx);
    
    
    for msIdx = 1:clusterMeasure_num
        clusterDescr(clIdx).([clusterMeasure_lbl{msIdx} '_threshold']) = clustThreshold.([clusterMeasure_lbl{msIdx} '_oneTail']);
        clusterDescr(clIdx).([clusterMeasure_lbl{msIdx} '_measure']) = clusterMeasure_obs.(clusterMeasure_lbl{msIdx})(clIdx);
        clusterDescr(clIdx).([clusterMeasure_lbl{msIdx} '_pvalue'])  = pval(clIdx,msIdx);
    end
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_0'  ]) = y1_0(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_25' ]) = y1_25(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_50' ]) = y1_50(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_75' ]) = y1_75(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_100']) = y1_100(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_range']) = ['[' num2str(y1_0(clIdx),'%.2f') ', ' num2str(y1_100(clIdx),'%.2f') ']'];
    
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_0'  ]) = x2_0(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_25' ]) = x2_25(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_50' ]) = x2_50(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_75' ]) = x2_75(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_100']) = x2_100(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_range']) = ['[' num2str(x2_0(clIdx),'%.2f') ', ' num2str(x2_100(clIdx),'%.2f') ']'];

    clusterDescr(clIdx).([DimStruct.z3_lbl '_count' ]) = chanInClusterNum(clIdx);
    clusterDescr(clIdx).([DimStruct.z3_lbl '_lbl' ])   = strjoin(chanInCluster{clIdx});

    clusterDescr(clIdx).eSiz_Mdn = eSizMdn(clIdx);

    clusterDescr(clIdx).eSiz_Ext = eSizExt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.y1_lbl '_' 'pnt' ])              = eSizExt_y1_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.y1_lbl '_' DimStruct.y1_units ]) = eSizExt_y1_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizExt_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizExt_x2_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.z3_lbl '_' 'lbl'])               = eSizExt_z3_lbl(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.z3_lbl '_' 'idx'])               = eSizExt_z3_idx(clIdx);

    clusterDescr(clIdx).eSiz_wMedoid = eSizMedoid(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.y1_lbl '_' 'pnt' ])              = eSizMedoid_y1_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.y1_lbl '_' DimStruct.y1_units ]) = eSizMedoid_y1_units(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizMedoid_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizMedoid_x2_units(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.z3_lbl '_' 'lbl'])               = eSizMedoid_z3_lbl(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.z3_lbl '_' 'idx'])               = eSizMedoid_z3_idx(clIdx);

    clusterDescr(clIdx).stat_Mdn = statMdn(clIdx);
    clusterDescr(clIdx).stat_Ext = statExt(clIdx);
end

% convert structure to table
clusterDescr = struct2table(clusterDescr);

% sort by criterion
sortBy = PE_parameters.clusterMetricChoice;  % TO DO: this is hard coded for now. give option to choose outside this script % DONE: THIS COMMENT CAN BE DELETED
if any(contains(clusterMeasure_lbl,sortBy))
    if nClust>0
        [~, sortIdx] = sort(abs(clusterDescr.([sortBy '_measure'])),'descend');
        clusterDescr = clusterDescr(sortIdx,:);
    end
else
    error('i cannot sort by the sortBy criterion: it does not exist')
end



