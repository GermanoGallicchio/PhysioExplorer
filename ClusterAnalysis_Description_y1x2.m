function clusterDescr = ...
    ClusterAnalysis_Description_y1x2(DimStruct, clusterMatrix_neg, clusterMatrix_pos, clusterMeasure_neg, clusterMeasure_pos, clustMaxMeasure_H0, statMatrix, eSizMatrix)

% ClusterAnalysis_Description_y1x2 describes the clusters 
%
% INPUT: 
%
% DimStruct         y1=Freq, x2=Time
%                   For example...
%                     DimStruct.y1_lbl   = 'Freq';
%                     DimStruct.y1_units = 'Hz';
%                     DimStruct.y1_vec   = CWT_f;
%                     DimStruct.x2_lbl   = 'Time';
%                     DimStruct.x2_units = 's';
%                     DimStruct.x2_vec   = activity_time;
%
% OUTPUT:
%
% 
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% debugging cell

% clustMaxMass_H0 = clustMaxMass1_H0;
%clusterMatrix_neg = clusterMatrix1_neg;
%clusterMatrix_pos = clusterMatrix1_pos;

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

ny1 = length(DimStruct.y1_vec);
nx2 = length(DimStruct.x2_vec);

[y1Matrix,       x2Matrix]       = ndgrid(1:ny1, 1:nx2);
[y1Matrix_units, x2Matrix_units] = ndgrid(DimStruct.y1_vec, DimStruct.x2_vec);

%%

% cluster measure (e.g., size, mass, rob mass)
clusterMeasure = [clusterMeasure_neg clusterMeasure_pos];
nClust = length(clusterMeasure);

% pvalues at cluster level (this corresponds with the inferential test)
% p value for each cluster (how many more extreme values are there in the H0 dist?)
% negative
pval_neg = nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    % BUG. DELETE THIS LINE pval_neg(clIdx) = sum(clustMaxMeasure_H0>abs(clusterMeasure_neg(clIdx))) / length(clustMaxMeasure_H0);
    pval_neg(clIdx) = sum(abs(clusterMeasure_neg(clIdx))<clustMaxMeasure_H0) / length(clustMaxMeasure_H0);
end
% positive
pval_pos = nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    % BUG. DELETE THIS LINE pval_pos(clIdx) = sum(clustMaxMeasure_H0>clusterMeasure_pos(clIdx)) / length(clustMaxMeasure_H0);
    pval_pos(clIdx) = sum(clusterMeasure_pos(clIdx)<clustMaxMeasure_H0) / length(clustMaxMeasure_H0);
end
% concatenate negative and positive
pval = [pval_neg pval_pos];


% y1 information
% find 0, 25, 50, 75, 100th prctile of the y1 dimension (e.g., freq) for each cluster
% negative
y1_0_neg  = nan(1,length(clusterMeasure_neg));
y1_25_neg = nan(1,length(clusterMeasure_neg));
y1_50_neg = nan(1,length(clusterMeasure_neg));
y1_75_neg = nan(1,length(clusterMeasure_neg));
y1_100_neg= nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    y1_0_neg(clIdx)  = prctile(y1Matrix_units(idx),0);
    y1_25_neg(clIdx) = prctile(y1Matrix_units(idx),25);
    y1_50_neg(clIdx) = median(y1Matrix_units(idx));
    y1_75_neg(clIdx) = prctile(y1Matrix_units(idx),75);
    y1_100_neg(clIdx)= prctile(y1Matrix_units(idx),100);
end
% positive
y1_0_pos  = nan(1,length(clusterMeasure_pos));
y1_25_pos = nan(1,length(clusterMeasure_pos));
y1_50_pos = nan(1,length(clusterMeasure_pos));
y1_75_pos = nan(1,length(clusterMeasure_pos));
y1_100_pos= nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    y1_0_pos(clIdx)  = prctile(y1Matrix_units(idx),0);
    y1_25_pos(clIdx) = prctile(y1Matrix_units(idx),25);
    y1_50_pos(clIdx) = median(y1Matrix_units(idx));
    y1_75_pos(clIdx) = prctile(y1Matrix_units(idx),75);
    y1_100_pos(clIdx)= prctile(y1Matrix_units(idx),100);
end
% concatenate negative and positive
y1_0  = [ y1_0_neg    y1_0_pos ];
y1_25 = [ y1_25_neg   y1_25_pos ];
y1_50 = [ y1_50_neg   y1_50_pos ];
y1_75 = [ y1_75_neg   y1_75_pos ];
y1_100= [ y1_100_neg  y1_100_pos ];


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





% effect size (e.g., Cohen's d, correlation coefficient)
% eSiz within each cluster (this does not correspond with any inferential test, purely descriptive)
%   - Mdn (median) 
%   - Ext (most extreme) of eSiz across all pixels within each cluster
%   - weighted Medoid
% negative
eSizMdn_neg = nan(1,length(clusterMeasure_neg));
eSizExt_neg = nan(1,length(clusterMeasure_neg));
eSizExt_y1_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizExt_y1_units_neg = nan(1,length(clusterMeasure_neg));
eSizExt_x2_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizExt_x2_units_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_y1_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_y1_units_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_x2_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_x2_units_neg = nan(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    % median
    eSizMdn_neg(clIdx) = median(eSizMatrix(idx));
    % most extreme
    [eSizExt_neg(clIdx),   ExtIdx] = min(eSizMatrix(idx));
    temp = y1Matrix(idx);           eSizExt_y1_pnt_neg(clIdx)   = temp(ExtIdx);
    temp = y1Matrix_units(idx);     eSizExt_y1_units_neg(clIdx) = temp(ExtIdx);
    temp = x2Matrix(idx);           eSizExt_x2_pnt_neg(clIdx)   = temp(ExtIdx);
    temp = x2Matrix_units(idx);     eSizExt_x2_units_neg(clIdx) = temp(ExtIdx);
    % weighted medoid 
    [~, ~, y1_wMedIdx, x2_wMedIdx] = ClusterAnalysis_y1x2Medoids(eSizMatrix,idx);
    eSizMedoid_y1_pnt_neg(clIdx)   = y1_wMedIdx;
    eSizMedoid_x2_pnt_neg(clIdx)   = x2_wMedIdx;
    eSizMedoid_y1_units_neg(clIdx) = y1Matrix_units(y1_wMedIdx,x2_wMedIdx);
    eSizMedoid_x2_units_neg(clIdx) = x2Matrix_units(y1_wMedIdx,x2_wMedIdx);
    eSizMedoid_neg(clIdx)  = eSizMatrix(y1_wMedIdx,x2_wMedIdx);
end
% positive
eSizMdn_pos = nan(1,length(clusterMeasure_pos));
eSizExt_pos = nan(1,length(clusterMeasure_pos));
eSizExt_y1_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizExt_y1_units_pos = nan(1,length(clusterMeasure_pos));
eSizExt_x2_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizExt_x2_units_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_y1_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_y1_units_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_x2_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_x2_units_pos = nan(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    % median
    eSizMdn_pos(clIdx) = median(eSizMatrix(idx));
    % most extreme
    [eSizExt_pos(clIdx),   ExtIdx] = max(eSizMatrix(idx));
    temp = y1Matrix(idx);           eSizExt_y1_pnt_pos(clIdx)   = temp(ExtIdx);
    temp = y1Matrix_units(idx);     eSizExt_y1_units_pos(clIdx) = temp(ExtIdx);
    temp = x2Matrix(idx);           eSizExt_x2_pnt_pos(clIdx)   = temp(ExtIdx);
    temp = x2Matrix_units(idx);     eSizExt_x2_units_pos(clIdx) = temp(ExtIdx);
    % weighted medoid 
    [~, ~, y1_wMedIdx, x2_wMedIdx] = ClusterAnalysis_y1x2Medoids(eSizMatrix,idx);
    eSizMedoid_y1_pnt_pos(clIdx)   = y1_wMedIdx;
    eSizMedoid_x2_pnt_pos(clIdx)   = x2_wMedIdx;
    eSizMedoid_y1_units_pos(clIdx) = y1Matrix_units(y1_wMedIdx,x2_wMedIdx);
    eSizMedoid_x2_units_pos(clIdx) = x2Matrix_units(y1_wMedIdx,x2_wMedIdx);
    eSizMedoid_pos(clIdx)  = eSizMatrix(y1_wMedIdx,x2_wMedIdx);
end
% concatenate negative and positive
eSizMdn = [ eSizMdn_neg   eSizMdn_pos ];
eSizExt = [ eSizExt_neg   eSizExt_pos ];
eSizExt_y1_pnt   = [ eSizExt_y1_pnt_neg    eSizExt_y1_pnt_pos   ];
eSizExt_y1_units = [ eSizExt_y1_units_neg  eSizExt_y1_units_pos ];
eSizExt_x2_pnt   = [ eSizExt_x2_pnt_neg    eSizExt_x2_pnt_pos   ];
eSizExt_x2_units = [ eSizExt_x2_units_neg  eSizExt_x2_units_pos ];
eSizMedoid = [ eSizMedoid_neg   eSizMedoid_pos ];
eSizMedoid_y1_pnt   = [ eSizMedoid_y1_pnt_neg    eSizMedoid_y1_pnt_pos  ];
eSizMedoid_y1_units = [ eSizMedoid_y1_units_neg  eSizMedoid_y1_units_pos];
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
    clusterDescr(clIdx).measure = clusterMeasure(clIdx);
    clusterDescr(clIdx).pvalue = pval(clIdx);

    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_0'  ]) = y1_0(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_25' ]) = y1_25(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_50' ]) = y1_50(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_75' ]) = y1_75(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_100']) = y1_100(clIdx);

    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_0'  ]) = x2_0(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_25' ]) = x2_25(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_50' ]) = x2_50(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_75' ]) = x2_75(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_100']) = x2_100(clIdx);

    clusterDescr(clIdx).eSiz_Mdn = eSizMdn(clIdx);

    clusterDescr(clIdx).eSiz_Ext = eSizExt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.y1_lbl '_' 'pnt' ])              = eSizExt_y1_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.y1_lbl '_' DimStruct.y1_units ]) = eSizExt_y1_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizExt_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Ext_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizExt_x2_units(clIdx);

    clusterDescr(clIdx).eSiz_wMedoid = eSizMedoid(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.y1_lbl '_' 'pnt' ])              = eSizMedoid_y1_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.y1_lbl '_' DimStruct.y1_units ]) = eSizMedoid_y1_units(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizMedoid_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_wMedoid_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizMedoid_x2_units(clIdx);

    clusterDescr(clIdx).stat_Mdn = statMdn(clIdx);
    clusterDescr(clIdx).stat_Ext = statExt(clIdx);
end

clusterDescr = struct2table(clusterDescr);
if nClust>0
    [~, sortIdx] = sort(abs(clusterDescr.measure),'descend');
    clusterDescr = clusterDescr(sortIdx,:);
end



