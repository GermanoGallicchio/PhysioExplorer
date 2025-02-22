function clusterDescr = ...
    ClusterAnalysis_Description_y1x2z3(DimStruct, clusterMatrix_neg, clusterMatrix_pos, clusterMeasure_neg, clusterMeasure_pos, clustMaxMeasure_H0, statMatrix, eSizMatrix)

% GG_oneDimPlusSpaceClusterDescription describes the clusters 
%
% INPUT: 
%
% DimStruct         y1=Freq, x2=Time, z=Chan
%                   For example...
%                     DimStruct.y1_lbl   = 'Freq';
%                     DimStruct.y1_units = 'Hz';
%                     DimStruct.y1_vec   = CWT_f;
%                     DimStruct.x2_lbl   = 'Time';
%                     DimStruct.x2_units = 's';
%                     DimStruct.x2_vec   = activity_time;
%                     DimStruct.z3_lbl   = 'Channel';
%                     DimStruct.z3_contFlag = 0;
%                     DimStruct.z3_chanlocs =; chanlocs structure in the eeglab format
%
% OUTPUT:
%
% 
% written by Germano Gallicchio 
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

if DimStruct.z3_contFlag==0
    %disp('z3 is not continuous (e.g., EEG channels)')
    chanlocs = DimStruct.z3_chanlocs;
else
    error('z3 is continuous. great but not yet done coded. if you only have two continuous variables (e.g., the other is channel), make sure that this is reflected in the DimStruct structure (i.e., the non-continuous variable should be z3)')
end

nFreq = length(DimStruct.y1_vec);
nTime = length(DimStruct.x2_vec);
nChan = length(chanlocs);

[z3Matrix,     y1Matrix,       x2Matrix]       = ndgrid(1:nChan,1:nFreq,1:nTime);
[z3Matrix_lbl, y1Matrix_units, x2Matrix_units] = ndgrid(string({chanlocs.labels}),DimStruct.y1_vec,DimStruct.x2_vec);

%%

% cluster measure (e.g., size, mass, rob mass)
clusterMeasure = [clusterMeasure_neg clusterMeasure_pos];
nClust = length(clusterMeasure);

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



% z3, channel information
% for each cluster, list channels representing all occurrences
% negative
chanInCluster_neg = cell(1,length(clusterMeasure_neg));
for clIdx = 1:length(clusterMeasure_neg)
    idx = clusterMatrix_neg==clIdx*-1;
    chanInCluster_neg{clIdx} = sort(unique(z3Matrix_lbl(idx)))';
end
% positive
chanInCluster_pos = cell(1,length(clusterMeasure_pos));
for clIdx = 1:length(clusterMeasure_pos)
    idx = clusterMatrix_pos==clIdx;
    chanInCluster_pos{clIdx} = sort(unique(z3Matrix_lbl(idx)))';
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
eSizExt_y1_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizExt_y1_units_neg = nan(1,length(clusterMeasure_neg));
eSizExt_x2_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizExt_x2_units_neg = nan(1,length(clusterMeasure_neg));
eSizExt_z3_lbl_neg   = repmat("",1,length(clusterMeasure_neg));
eSizExt_z3_idx_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_y1_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_y1_units_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_x2_pnt_neg   = nan(1,length(clusterMeasure_neg));
eSizMedoid_x2_units_neg = nan(1,length(clusterMeasure_neg));
eSizMedoid_z3_lbl_neg   = repmat("",1,length(clusterMeasure_neg));
eSizMedoid_z3_idx_neg   = nan(1,length(clusterMeasure_neg));
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
    temp = z3Matrix_lbl(idx);       eSizExt_z3_lbl_neg(clIdx)   = temp(ExtIdx);
    temp = z3Matrix(idx);           eSizExt_z3_idx_neg(clIdx)   = temp(ExtIdx);
    % weighted medoid 
    % compute dissimilarity between each pairs in the cluster (points a and b) as Euclidean distance, using binned angular distance for z3
    eSiz_clustPointsVals = eSizMatrix(idx);
    y1_clustPointsCoords = y1Matrix(idx); % y1-dimension coordinates for all points in the cluster
    x2_clustPointsCoords = x2Matrix(idx); % x2-dimension "
    z3_clustPointsCoords = z3Matrix(idx); % z3-dimension "
    y1_clustPointsUnits = y1Matrix_units(idx); % y1-dimension units for all points in the cluster
    x2_clustPointsUnits = x2Matrix_units(idx); % x2-dimension "
    z3_clustPointsLbls  = z3Matrix_lbl(idx);   % z3-dimension labels "
    euclDistance = nan(nnz(idx),nnz(idx));
    for aIdx = 1:nnz(idx) 
        for bIdx = 1:nnz(idx)
            euclDistance(aIdx,bIdx) = sqrt(  (y1_clustPointsCoords(aIdx)-y1_clustPointsCoords(bIdx)).^2  +  (x2_clustPointsCoords(aIdx)-x2_clustPointsCoords(bIdx)).^2  +  (DimStruct.z3_Dbinned(z3Matrix(aIdx),z3Matrix(bIdx))).^2  ); 
        end
    end % TO DO: vectorize the following lines to avoid for loops and achieve faster computation time
    [~, MedoidIdx] = min(sum(euclDistance,1) .* abs(eSiz_clustPointsVals').^-1);
    eSizMedoid_y1_pnt_neg(clIdx)   = y1_clustPointsCoords(MedoidIdx);
    eSizMedoid_y1_units_neg(clIdx) = y1_clustPointsUnits(MedoidIdx);
    eSizMedoid_x2_pnt_neg(clIdx)   = x2_clustPointsCoords(MedoidIdx);
    eSizMedoid_x2_units_neg(clIdx) = x2_clustPointsUnits(MedoidIdx);
    eSizMedoid_z3_lbl_neg(clIdx)   = z3_clustPointsLbls(MedoidIdx);
    eSizMedoid_z3_idx_neg(clIdx)   = z3_clustPointsCoords(MedoidIdx);
    eSizMedoid_neg(clIdx)  = eSiz_clustPointsVals(MedoidIdx);
end
% positive
eSizMdn_pos = nan(1,length(clusterMeasure_pos));
eSizExt_pos = nan(1,length(clusterMeasure_pos));
eSizExt_y1_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizExt_y1_units_pos = nan(1,length(clusterMeasure_pos));
eSizExt_x2_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizExt_x2_units_pos = nan(1,length(clusterMeasure_pos));
eSizExt_z3_lbl_pos   = repmat("",1,length(clusterMeasure_pos));
eSizExt_z3_idx_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_y1_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_y1_units_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_x2_pnt_pos   = nan(1,length(clusterMeasure_pos));
eSizMedoid_x2_units_pos = nan(1,length(clusterMeasure_pos));
eSizMedoid_z3_lbl_pos   = repmat("",1,length(clusterMeasure_pos));
eSizMedoid_z3_idx_pos   = nan(1,length(clusterMeasure_pos));
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
    temp = z3Matrix_lbl(idx);       eSizExt_z3_lbl_pos(clIdx)   = temp(ExtIdx);
    temp = z3Matrix(idx);           eSizExt_z3_idx_pos(clIdx)   = temp(ExtIdx);
    % weighted medoid 
    % compute dissimilarity between each pairs in the cluster (points a and b) as Euclidean distance, using binned angular distance for z3
    eSiz_clustPointsVals = eSizMatrix(idx);
    y1_clustPointsCoords = y1Matrix(idx); % y1-dimension coordinates for all points in the cluster
    x2_clustPointsCoords = x2Matrix(idx); % x2-dimension "
    z3_clustPointsCoords = z3Matrix(idx); % z3-dimension "
    y1_clustPointsUnits = y1Matrix_units(idx); % y1-dimension units for all points in the cluster
    x2_clustPointsUnits = x2Matrix_units(idx); % x2-dimension "
    z3_clustPointsLbls  = z3Matrix_lbl(idx);   % z3-dimension labels "
    euclDistance = nan(nnz(idx),nnz(idx));
    for aIdx = 1:nnz(idx) 
        for bIdx = 1:nnz(idx)
            euclDistance(aIdx,bIdx) = sqrt(  (y1_clustPointsCoords(aIdx)-y1_clustPointsCoords(bIdx)).^2  +  (x2_clustPointsCoords(aIdx)-x2_clustPointsCoords(bIdx)).^2  +  (DimStruct.z3_Dbinned(z3Matrix(aIdx),z3Matrix(bIdx))).^2  ); 
        end
    end % TO DO: vectorize the following lines to avoid for loops and achieve faster computation time
    [~, MedoidIdx] = min(sum(euclDistance,1) .* eSiz_clustPointsVals'.^-1);
    eSizMedoid_y1_pnt_pos(clIdx)   = y1_clustPointsCoords(MedoidIdx);
    eSizMedoid_y1_units_pos(clIdx) = y1_clustPointsUnits(MedoidIdx);
    eSizMedoid_x2_pnt_pos(clIdx)   = x2_clustPointsCoords(MedoidIdx);
    eSizMedoid_x2_units_pos(clIdx) = x2_clustPointsUnits(MedoidIdx);
    eSizMedoid_z3_lbl_pos(clIdx)   = z3_clustPointsLbls(MedoidIdx);
    eSizMedoid_z3_idx_pos(clIdx)   = z3_clustPointsCoords(MedoidIdx);
    eSizMedoid_pos(clIdx)  = eSiz_clustPointsVals(MedoidIdx);
end
% concatenate negative and positive
eSizMdn = [ eSizMdn_neg   eSizMdn_pos ];
eSizExt = [ eSizExt_neg   eSizExt_pos ];
eSizExt_y1_pnt   = [ eSizExt_y1_pnt_neg    eSizExt_y1_pnt_pos   ];
eSizExt_y1_units = [ eSizExt_y1_units_neg  eSizExt_y1_units_pos ];
eSizExt_x2_pnt   = [ eSizExt_x2_pnt_neg    eSizExt_x2_pnt_pos   ];
eSizExt_x2_units = [ eSizExt_x2_units_neg  eSizExt_x2_units_pos ];
eSizExt_z3_lbl   = [ eSizExt_z3_lbl_neg    eSizExt_z3_lbl_pos   ];
eSizExt_z3_idx   = [ eSizExt_z3_idx_neg    eSizExt_z3_idx_pos   ];
eSizMedoid = [ eSizMedoid_neg   eSizMedoid_pos ];
eSizMedoid_y1_pnt   = [ eSizMedoid_y1_pnt_neg    eSizMedoid_y1_pnt_pos  ];
eSizMedoid_y1_units = [ eSizMedoid_y1_units_neg  eSizMedoid_y1_units_pos];
eSizMedoid_x2_pnt   = [ eSizMedoid_x2_pnt_neg    eSizMedoid_x2_pnt_pos  ];
eSizMedoid_x2_units = [ eSizMedoid_x2_units_neg  eSizMedoid_x2_units_pos];
eSizMedoid_z3_lbl   = [ eSizMedoid_z3_lbl_neg    eSizMedoid_z3_lbl_pos  ];
eSizMedoid_z3_idx   = [ eSizMedoid_z3_idx_neg    eSizMedoid_z3_idx_pos  ];




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

    clusterDescr(clIdx).([DimStruct.y1_lbl 'y1_' DimStruct.y1_units '_0'  ]) = y1_0(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_25' ]) = y1_25(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_50' ]) = y1_50(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_75' ]) = y1_75(clIdx);
    clusterDescr(clIdx).([DimStruct.y1_lbl '_' DimStruct.y1_units '_100']) = y1_100(clIdx);

    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_0'  ]) = x2_0(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_25' ]) = x2_25(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_50' ]) = x2_50(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_75' ]) = x2_75(clIdx);
    clusterDescr(clIdx).([DimStruct.x2_lbl '_' DimStruct.x2_units '_100']) = x2_100(clIdx);

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

    clusterDescr(clIdx).eSiz_Medoid = eSizMedoid(clIdx);
    clusterDescr(clIdx).(['eSiz_Medoid_' DimStruct.y1_lbl '_' 'pnt' ])              = eSizMedoid_y1_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Medoid_' DimStruct.y1_lbl '_' DimStruct.y1_units ]) = eSizMedoid_y1_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Medoid_' DimStruct.x2_lbl '_' 'pnt' ])              = eSizMedoid_x2_pnt(clIdx);
    clusterDescr(clIdx).(['eSiz_Medoid_' DimStruct.x2_lbl '_' DimStruct.x2_units ]) = eSizMedoid_x2_units(clIdx);
    clusterDescr(clIdx).(['eSiz_Medoid_' DimStruct.z3_lbl '_' 'lbl'])               = eSizMedoid_z3_lbl(clIdx);
    clusterDescr(clIdx).(['eSiz_Medoid_' DimStruct.z3_lbl '_' 'idx'])               = eSizMedoid_z3_idx(clIdx);

    clusterDescr(clIdx).stat_Mdn = statMdn(clIdx);
    clusterDescr(clIdx).stat_Ext = statExt(clIdx);
end

clusterDescr = struct2table(clusterDescr);
if nClust>0
    [~, sortIdx] = sort(abs(clusterDescr.measure),'descend');
    clusterDescr = clusterDescr(sortIdx,:);
end



