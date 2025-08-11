function results = pe_describeClusters(pe_cfg,results)

% Descr
%
% permutation (hypothesis-testing inference) => pvalues (1 tail)
% bootstrap (estimation-based inference) => BR, 95%CI (2 tails)
%
% permutation gives a 1-tail significanceThreshold
% bootstrap gives a 2-tail range 
%
% INPUT:
%
%   pe_cfg
%
%
% OUTPUT:
%
% rowIdx
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

[y1Matrix_units, x2Matrix_units, z3Matrix_lbl ] = ndgrid(pe_cfg.dimensions.y1_vec,  pe_cfg.dimensions.x2_vec,  pe_cfg.dimensions.z3_chanLbl);

y1Lbl = pe_cfg.dimensions.y1_lbl;
y1Units = pe_cfg.dimensions.y1_units;
x2Lbl = pe_cfg.dimensions.x2_lbl;
x2Units = pe_cfg.dimensions.x2_units;
z3Lbl = pe_cfg.dimensions.z3_lbl;

%% descriptive clusters 
% where no inferential cluster analysis happened already
% (e.g., any of the bootstrap analyses, empiricalL1_FDR analysis)

switch [pe_cfg.objective ' & ' pe_cfg.analysis]
    case 'permutationH0testing & empiricalL1_FDR'
        % identify clusters on empirical pvalue results
        mask = results.resampling.pVal_emp_FDR<pe_cfg.p_crit;  % binarize the pvalue
        mask = mask .* ~pe_cfg.R_ignore;  % apply R_ignore mask

        % cluster forming
        [clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, results.statVal_obs , ~mask);
        results.clusters.clusterMembership_obs = clusterMembership;
        results.clusters.clustIDList_obs       = clustIDList;
        results.clusters.metrics_obs           = metrics;

    case 'bootstrapStability & empiricalL1_FDR'
        % binarize the bootstrap results and find and describe clusters
        mask_BRrob = logical(abs(results.resampling.inference.BR_rob) > 2);
        mask_CI    = logical(abs(sum(sign([results.resampling.inference.CIlo; results.resampling.inference.CIup]),1))>0);
        % apply R_ignore mask
        mask_BRrob = mask_BRrob .* ~pe_cfg.R_ignore;
        mask_CI    = mask_CI    .* ~pe_cfg.R_ignore;

        % clusters based on bootstrap criteria (to choose below by muting either line)
        [clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, results.statVal_obs , ~mask_BRrob);
        %[clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, results.statVal_obs , ~mask_CI);
        results.clusters.clusterMembership_obs = clusterMembership;
        results.clusters.clustIDList_obs       = clustIDList;
        results.clusters.metrics_obs           = metrics;

    case 'permutationH0testing & theoreticalL1_clusterMaxT'
        % nothing to do / clusters were already created as main idea of this analysis

    case 'bootstrapStability & theoreticalL1_clusterMaxT'
        % identify clusters on results of bootstrap
        warning('double-check this works'); keyboard

        % binarize the bootstrap results and find and describe clusters
        mask_BRrob = logical(abs(results.lvl1.inference_maxT.BR_rob) > 2);
        mask_CI    = logical(abs(sum(sign([results.lvl1.inference_maxT.CIlo; results.lvl1.inference_maxT.CIup]),1))/2);

        % apply R_ignore mask
        mask_BRrob = mask_BRrob .* ~pe_cfg.R_ignore;
        mask_CI    = mask_CI    .* ~pe_cfg.R_ignore;

        % clusters based on bootstrap criteria (to choose below by muting either line)
        [clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, results.lvl1.statVal , ~mask_BRrob);
        %[clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, results.lvl1.statVal , ~mask_CI);
        results.clusters.clusterMembership = clusterMembership;
        results.clusters.clustIDList       = clustIDList;
        results.clusters.metrics           = metrics;

    case 'permutationH0testing & PLS_SVD'
        % nothing to do / no clusters to describe
    
    case 'bootstrapStability & PLS_SVD'
        error('not coded yet')
end

%% more shortcuts

clusterMembership = results.clusters.clusterMembership_obs;
clustIDList = results.clusters.clustIDList_obs;
nClust = length(clustIDList);

%% cluster descriptives

% initialize vars
statVal_Mdn        = nan(1,nClust);
statVal_IQR        = nan(1,nClust);
statVal_mostExtrm  = nan(1,nClust);
y1_units_mostExtrm = nan(1,nClust);
x2_units_mostExtrm = nan(1,nClust);
z3_lbl_mostExtrm   = repmat("",1,nClust);
statVal_leastExtrm  = nan(1,nClust);
y1_units_leastExtrm = nan(1,nClust);
x2_units_leastExtrm = nan(1,nClust);
z3_lbl_leastExtrm   = repmat("",1,nClust);
inferVal_Mdn        = nan(1,nClust);
inferVal_IQR        = nan(1,nClust);
inferVal_mostExtrm  = nan(1,nClust);
inferVal_leastExtrm = nan(1,nClust);
y1_num       = nan(1,nClust);
y1_units_min = nan(1,nClust);
y1_units_p25 = nan(1,nClust);
y1_units_p75 = nan(1,nClust);
y1_units_max = nan(1,nClust);
x2_num       = nan(1,nClust);
x2_units_min = nan(1,nClust);
x2_units_p25 = nan(1,nClust);
x2_units_p75 = nan(1,nClust);
x2_units_max = nan(1,nClust);
z3_num       = nan(1,nClust);
z3_lbl       = repmat("",1,nClust);
for clIdx = 1:nClust
    
    % find idx corresponding with this cluster
    idx = clusterMembership==clustIDList(clIdx);
    
    % num, range, and iqr (metrics of extension) along each dimension
    y1Vals = unique(y1Matrix_units(idx));
    y1Vals(isnan(y1Vals)) = []; y1Vals = [y1Vals NaN]; % in case of multiple NaNs, just keep one NaN
    y1_num(1,clIdx)       = length(y1Vals);     % number of unique values
    y1_units_min(1,clIdx) = min(y1Vals);        % min among unique values
    y1_units_p25(1,clIdx) = prctile(y1Vals,25); % 25th prctle among unique values
    y1_units_p75(1,clIdx) = prctile(y1Vals,75); % 75th prctle among unique values
    y1_units_max(1,clIdx) = max(y1Vals);        % max among unique values
    x2Vals = unique(x2Matrix_units(idx));
    x2Vals(isnan(x2Vals)) = []; x2Vals = [x2Vals NaN]; % in case of multiple NaNs, just keep one NaN
    x2_num(1,clIdx)       = length(x2Vals);
    x2_units_min(1,clIdx) = min(x2Vals);
    x2_units_p25(1,clIdx) = prctile(x2Vals,25);
    x2_units_p75(1,clIdx) = prctile(x2Vals,75);
    x2_units_max(1,clIdx) = max(x2Vals);
    z3_num(1,clIdx) = length(unique(z3Matrix_lbl(idx)));
    z3_lbl(1,clIdx) = join(unique(z3Matrix_lbl(idx)));

    % all statVals within this cluster
    vals = results.statVal_obs(idx);

    % median statVal
    statVal_Mdn(1,clIdx) = median(vals);

    % iqr statVal
    statVal_IQR(1,clIdx) = iqr(vals);

    % most extreme statVal
    [~, mostExtrm_idx] = max(abs(vals));
    statVal_mostExtrm(1,clIdx) = vals(mostExtrm_idx);

    % most extreme location
    temp = y1Matrix_units(idx);     
    y1_units_mostExtrm(1,clIdx) = temp(mostExtrm_idx);
    temp = x2Matrix_units(idx);     
    x2_units_mostExtrm(1,clIdx) = temp(mostExtrm_idx);
    temp = z3Matrix_lbl(idx);     
    z3_lbl_mostExtrm(1,clIdx) = temp(mostExtrm_idx);
    
    % least extreme statVal
    [~, leastExtrm_idx] = min(abs(vals));
    statVal_leastExtrm(1,clIdx) = vals(leastExtrm_idx);

    % least extreme location
    temp = y1Matrix_units(idx);     
    y1_units_leastExtrm(1,clIdx) = temp(leastExtrm_idx);
    temp = x2Matrix_units(idx);     
    x2_units_leastExtrm(1,clIdx) = temp(leastExtrm_idx);
    temp = z3Matrix_lbl(idx);     
    z3_lbl_leastExtrm(1,clIdx) = temp(leastExtrm_idx);
    
    % bootstrap ratio values within this cluster
    switch pe_cfg.objective
        case 'bootstrapStability'
            vals = results.resampling.inference.BR_rob(idx);

            % median inference metric
            inferVal_Mdn(1,clIdx) = median(vals);

            % iqr inference metric
            inferVal_IQR(1,clIdx) = iqr(vals);

            % most extreme inference metric
            [~, mostExtrm_idx] = max(abs(vals));
            inferVal_mostExtrm(1,clIdx) = vals(mostExtrm_idx);

            % least extreme inference metric
            [~, leastExtrm_idx] = min(abs(vals));
            inferVal_leastExtrm(1,clIdx) = vals(leastExtrm_idx);
    end

    
    
end
results.clusters.descriptive.statVal_Mdn = statVal_Mdn;
results.clusters.descriptive.statVal_IQR = statVal_IQR;

results.clusters.descriptive.statVal_mostExtrm = statVal_mostExtrm;
results.clusters.descriptive.([y1Lbl '_' y1Units '_mostExtrm']) = y1_units_mostExtrm;
results.clusters.descriptive.([x2Lbl '_' x2Units '_mostExtrm']) = x2_units_mostExtrm;
results.clusters.descriptive.([z3Lbl '_mostExtrm'])             = z3_lbl_mostExtrm;

results.clusters.descriptive.statVal_leastExtrm = statVal_leastExtrm;
results.clusters.descriptive.([y1Lbl '_' y1Units '_leastExtrm']) = y1_units_leastExtrm;
results.clusters.descriptive.([x2Lbl '_' x2Units '_leastExtrm']) = x2_units_leastExtrm;
results.clusters.descriptive.([z3Lbl '_leastExtrm'])             = z3_lbl_leastExtrm;

if strcmp(pe_cfg.objective,'bootstrapStability')
    results.clusters.descriptive.(['BRrob_Mdn'])        = inferVal_Mdn;
    results.clusters.descriptive.(['BRrob_IQR'])        = inferVal_IQR;
    results.clusters.descriptive.(['BRrob_mostExtrm'])  = inferVal_mostExtrm;
    results.clusters.descriptive.(['BRrob_leastExtrm']) = inferVal_leastExtrm;
end
results.clusters.descriptive.([y1Lbl '_num' ])       = y1_num;                
results.clusters.descriptive.([y1Lbl '_' y1Units '_p00']) = y1_units_min;
results.clusters.descriptive.([y1Lbl '_' y1Units '_p25']) = y1_units_p25;
results.clusters.descriptive.([y1Lbl '_' y1Units '_p75']) = y1_units_p75;
results.clusters.descriptive.([y1Lbl '_' y1Units '_p100']) = y1_units_max;
results.clusters.descriptive.([x2Lbl '_num' ])       = x2_num;
results.clusters.descriptive.([x2Lbl '_' x2Units '_p00']) = x2_units_min;
results.clusters.descriptive.([x2Lbl '_' x2Units '_p25']) = x2_units_p25;
results.clusters.descriptive.([x2Lbl '_' x2Units '_p75']) = x2_units_p75;
results.clusters.descriptive.([x2Lbl '_' x2Units '_p100']) = x2_units_max;
results.clusters.descriptive.([z3Lbl '_num' ])       = z3_num;
results.clusters.descriptive.([z3Lbl '_lbl' ])       = z3_lbl;

%% summary table
summaryStruct = struct(); % initialize
for clIdx = 1:nClust
    % metrics
    tmp = results.clusters.metrics_obs;
    fieldLbl = fieldnames(tmp);
    fieldNum = length(fieldLbl);
    for fIdx = 1:fieldNum
        summaryStruct(clIdx).(fieldLbl{fIdx}) = tmp.(fieldLbl{fIdx})(clIdx);
    end
    % inference / pvalues
    if strcmp(pe_cfg.objective,'permutationH0testing')  &  strcmp(pe_cfg.analysis,'theoreticalL1_clusterMaxT')
        tmp = results.clusters.inference_maxT.pval;
        fieldLbl = fieldnames(tmp);
        fieldNum = length(fieldLbl);
        for fIdx = 1:fieldNum
            summaryStruct(clIdx).(['pval_' fieldLbl{fIdx}]) = tmp.(fieldLbl{fIdx})(clIdx);
        end
    end
    % inference / significance thresholds
    if contains(pe_cfg.objective,'permutationH0testing')  &  strcmp(pe_cfg.analysis,'theoreticalL1_clusterMaxT')
        tmp = results.clusters.inference_maxT.thresholds;
        fieldLbl = fieldnames(tmp);
        fieldNum = length(fieldLbl);
        for fIdx = 1:fieldNum
            summaryStruct(clIdx).(['thresholds_' fieldLbl{fIdx}]) = tmp.(fieldLbl{fIdx});
        end
    end
    % descriptive
    tmp = results.clusters.descriptive;
    fieldLbl = fieldnames(tmp);
    fieldNum = length(fieldLbl);
    for fIdx = 1:fieldNum
        summaryStruct(clIdx).(fieldLbl{fIdx}) = tmp.(fieldLbl{fIdx})(clIdx);
    end
end
summaryTable = struct2table(summaryStruct);

% sort rows by cluster metrics combined
metrics_lbl = fieldnames(results.clusters.metrics_obs);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id')); % remove the 'id' column
metrics_num = length(metrics_lbl);
measureVec = zeros(1,nClust);
for msIdx = 1:metrics_num
    measureVec = measureVec + normalize(results.clusters.metrics_obs.(metrics_lbl{msIdx}),'range',[0 1]);
end
if nClust>0
    [~, sortIdx] = sort(measureVec,'descend');
    summaryTable = summaryTable(sortIdx,:);
end
results.clusters.descriptive.summaryTable = summaryTable;

% if pe_cfg.verbose   
%     disp('level-2 descriptives added')
% end
