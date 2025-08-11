function results = pe_description(pe_cfg, results)

% PE_Description_y1x2z3 describes the results
%
% INPUT: 
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
%                        - All fields are labeled using the provided pe_cfg.dimensions labels/units
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

%% shortcuts

ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;



statVal           = results.featureWise.statVal;
clustIDList       = results.featureWise.clustIDList;
clusterMembership = results.featureWise.clusterMembership;
metrics           = results.clusterWise.metrics; 
signifThreshold   = results.clusterWise.signifThreshold;
statVal_y1x2z3 = reshape(statVal,[ny1 nx2 nz3]); % convert to y1-x2-z3 shape
clusterMembership_y1x2z3 = reshape(clusterMembership,[ny1 nx2 nz3]); % convert to y1-x2-z3 shape


nClust = length(clustIDList);


metrics_lbl = fieldnames(metrics);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id')); % remove the 'id' column
metrics_num = length(metrics_lbl);

%%

[y1Matrix,       x2Matrix,       z3Matrix,    ] = ndgrid(1:ny1,1:nx2,1:nz3);
[y1Matrix_units, x2Matrix_units, z3Matrix_lbl ] = ndgrid(pe_cfg.dimensions.y1_vec,pe_cfg.dimensions.x2_vec,string({pe_cfg.dimensions.z3_chanLbl}));




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
    idx = clusterMembership==clustIDList(clIdx);
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
    idx = clusterMembership==clustIDList(clIdx);
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
    idx = clusterMembership==clustIDList(clIdx);
    chanInCluster{clIdx} = unique(z3Matrix_lbl(idx))';
end
chanInClusterNum = cellfun(@length,chanInCluster);



% effect size (e.g., Cohen's d, correlation coefficient)
% eSiz within each cluster (this does not correspond with any inferential test, purely descriptive)
%   - Mdn (median) 
%   - Ext (most extreme) of eSiz across all pixels within each cluster
eSizMdn = nan(1,nClust);
eSizExt = nan(1,nClust);
eSizExt_y1_pnt   = nan(1,nClust);
eSizExt_y1_units = nan(1,nClust);
eSizExt_x2_pnt   = nan(1,nClust);
eSizExt_x2_units = nan(1,nClust);
eSizExt_z3_lbl   = repmat("",1,nClust);
eSizExt_z3_idx   = nan(1,nClust);
for clIdx = 1:nClust
    idx = clusterMembership==clustIDList(clIdx);
    % median
    eSizMdn(clIdx) = median(statVal(idx));
    % most extreme
    if sign(clustIDList(clIdx)) < 0
        [eSizExt(clIdx),   ExtIdx] = min(statVal(idx)); % min for negative clusters
    else
        [eSizExt(clIdx),   ExtIdx] = max(statVal(idx)); % max for positive clusters
    end
    temp = y1Matrix(idx);           eSizExt_y1_pnt(clIdx)   = temp(ExtIdx);
    temp = y1Matrix_units(idx);     eSizExt_y1_units(clIdx) = temp(ExtIdx);
    temp = x2Matrix(idx);           eSizExt_x2_pnt(clIdx)   = temp(ExtIdx);
    temp = x2Matrix_units(idx);     eSizExt_x2_units(clIdx) = temp(ExtIdx);
    temp = z3Matrix_lbl(idx);       eSizExt_z3_lbl(clIdx)   = temp(ExtIdx);
    temp = z3Matrix(idx);           eSizExt_z3_idx(clIdx)   = temp(ExtIdx);
end


% stat values at cluster level (this does not correspond with any test, purely descriptive)
% Mdn (median) and Ext (most extreme) of stat (e.g., tvalue, corr coefficients) across all pixels within each cluster
statMdn = nan(1,nClust);
statExt = nan(1,nClust);
for clIdx = 1:nClust
    idx = clusterMembership==clustIDList(clIdx);
    statMdn(clIdx) = median(statVal(idx));
    if sign(clustIDList(clIdx)) < 0
        statExt(clIdx) = min(statVal(idx));  % min for negative
    else
        statExt(clIdx) = max(statVal(idx));  % max for positive
    end
    
    
end



clusterDescr = struct();
for clIdx = 1:nClust
    clusterDescr(clIdx).id        = clustIDList(clIdx);
    
    
    for msIdx = 1:metrics_num
        clusterDescr(clIdx).([metrics_lbl{msIdx} '_threshold']) = signifThreshold.([metrics_lbl{msIdx}]);
        clusterDescr(clIdx).([metrics_lbl{msIdx} '_measure']) = metrics.(metrics_lbl{msIdx})(clIdx);
        clusterDescr(clIdx).([metrics_lbl{msIdx} '_pvalue'])  = results.clusterWise.inference.(['pval_' metrics_lbl{msIdx}])(clIdx);
    end
    clusterDescr(clIdx).([pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units '_0'  ]) = y1_0(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units '_25' ]) = y1_25(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units '_50' ]) = y1_50(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units '_75' ]) = y1_75(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units '_100']) = y1_100(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units '_range']) = ['[' num2str(y1_0(clIdx),'%.2f') ', ' num2str(y1_100(clIdx),'%.2f') ']'];
    
    clusterDescr(clIdx).([pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units '_0'  ]) = x2_0(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units '_25' ]) = x2_25(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units '_50' ]) = x2_50(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units '_75' ]) = x2_75(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units '_100']) = x2_100(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units '_range']) = ['[' num2str(x2_0(clIdx),'%.2f') ', ' num2str(x2_100(clIdx),'%.2f') ']'];

    clusterDescr(clIdx).([pe_cfg.dimensions.z3_lbl '_count' ]) = chanInClusterNum(clIdx);
    clusterDescr(clIdx).([pe_cfg.dimensions.z3_lbl '_lbl' ])   = strjoin(chanInCluster{clIdx});

    clusterDescr(clIdx).eSiz_Mdn = eSizMdn(clIdx);

    clusterDescr(clIdx).eSiz_Ext = eSizExt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' pe_cfg.dimensions.y1_lbl '_' 'pnt' ])              = eSizExt_y1_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' pe_cfg.dimensions.y1_lbl '_' pe_cfg.dimensions.y1_units ]) = eSizExt_y1_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' pe_cfg.dimensions.x2_lbl '_' 'pnt' ])              = eSizExt_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' pe_cfg.dimensions.x2_lbl '_' pe_cfg.dimensions.x2_units ]) = eSizExt_x2_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' pe_cfg.dimensions.z3_lbl '_' 'lbl'])               = eSizExt_z3_lbl(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' pe_cfg.dimensions.z3_lbl '_' 'idx'])               = eSizExt_z3_idx(clIdx);

    clusterDescr(clIdx).stat_Mdn = statMdn(clIdx);
    clusterDescr(clIdx).stat_Ext = statExt(clIdx);
end

% convert structure to table
clusterDescr = struct2table(clusterDescr);

% sort by criterion
% sort by average of p values across all measures

if nClust>0
    [~, sortIdx] = sort(mean(cell2mat(struct2cell(results.clusterWise.inference)),1));
    clusterDescr = clusterDescr(sortIdx,:);
end

results.clusterWise.clusterDescr = clusterDescr;

