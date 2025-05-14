function [results, clusterMetrics, clustThreshold] = ...
    PE_Stats_y1x2z3(dataArray, ...
    PE_parameters, ...
    DimStruct)
% Performs statistical analysis on EEG data with cluster-based inference
%
% [results, clusterMetrics, clustThreshold] = PE_Stats_y1x2z3(dataArray, PE_parameters, DimStruct)
%
% This function performs statistical analysis on EEG data using cluster-based inference methods.
% It supports group comparisons, condition comparisons, and correlation analyses with options
% for different statistical tests and inference methods.

% INPUT:
%   dataArray     - 3D array of data: subjects × dimensions × (optional) conditions
%
%   PE_parameters - Structure containing analysis parameters:
%     .test         - Type of test: 'group2vs1_ttest', 'condition2vs1_ttest',
%                     'correlation_Pearson', 'correlation_Spearman', or 'correlation_Kendall'
%     .stats        - 'H0MonteCarlo' for permutation/bootstrap or 'observed' for single analysis
%     .inference    - 'Tmax' for maximum statistic (recommended) or 'Tall' for all clusters (not recommended)
%     .cluster_p_crit - Critical p-value for clustering
%     .group        - Vector of group assignments for subjects (all 1s if there is only one group--if only interested in a subset of the sample, change dataArray accordingly)
%     .H0MonteCarlo_replacement - Boolean for bootstrap (true) vs permutation (false)
%     .nIterations  - Number of iterations for Monte Carlo methods
%     .randomSeed   - Seed for random number generator
%     .verbose      - Boolean for verbose output
%
%   DimStruct    - Structure defining data dimensions / same as in earlier functions (e.g., PE_Adjacency_y1x2z3)
%
%
% OUTPUT:
%
%   results        - Structure containing test results (depending on the statistical test):
%     .tvals        - t-values (for t-tests)
%     .dCohen       - Cohen's d effect sizes (for t-tests)
%     .rvals        - Correlation coefficients (for correlation analyses)
%     .pvals        - p-values
%     .clusterMatrix - Matrix identifying clusters
%     .clustIDList  - List of cluster IDs
%
%   clusterMetrics - Structure array with cluster statistics:
%     .id          - Cluster identifiers
%     .size        - Cluster sizes
%     .mass        - Cluster masses
%
%   clustThreshold - Structure containing threshold values for cluster inference:
%     .size_oneTail   - One-tailed threshold for cluster size
%     .size_twoTails  - Two-tailed thresholds for cluster size
%     .mass_oneTail   - One-tailed threshold for cluster mass
%     .mass_twoTails  - Two-tailed thresholds for cluster mass
%
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

dataArray_orig  = dataArray;

if strcmp(PE_parameters.inference,'Tall')
    if strcmp(PE_parameters.stats, 'H0MonteCarlo')
        warning('Tall is exploratory. It does not control for FWER whereas Tmax does')
    end
end

% PE_parameters.stats label
switch PE_parameters.stats
    case {'H0MonteCarlo'  'observed'}
        if PE_parameters.verbose
            disp(['using: ' PE_parameters.stats])
        end
    otherwise
        error('choose either H0MonteCarlo or observed')
end



% number of iterations
switch PE_parameters.stats
    case 'H0MonteCarlo'
        nIterations = PE_parameters.nIterations;
    case 'observed'
        nIterations = 1; % no iterations
end

%% parameters
% TO DO: sanity check that participants are sorted by group information

nSubj = size(dataArray,1);

% group comparison
if strcmp(PE_parameters.test,'group2vs1_ttest')
    if any(contains(fieldnames(PE_parameters),'group'))
        group = PE_parameters.group;
        groupLbl = unique(group);
        groupNum = length(groupLbl);

        % numerosity of each group
        countTbl = countlabels(group);
        countTbl = sortrows(countTbl,"Label",'ascend'); % enforce ascending order
        if ~isequal(double(string(countTbl.Label)),groupLbl); warning('error with grouping variable'); keyboard; end

    else
        error('need to specify the PE_parameters.group vector for this analysis')
    end
    if ~isequal(length(group), size(dataArray,1))
        error('the group variable needs to be same length as the size(dataArray,1)')
    end
    if any(contains(fieldnames(PE_parameters),'condition'))
        condition = PE_parameters.condition;
        if size(dataArray,3)~=length(condition)
            error('condition vector and size(dataArray,3) should match')
        end
        if ~any(condition==1)
            error('include at least one ''1'' in the condition vector')
            if any(condition~=1  &  condition~=0)
                error('include ''1'' and ''0'' in the condition vector to label those to be included')
            end
        end
        nConditions = sum(condition==1);
    end
end

% condition comparison
if strcmp(PE_parameters.test,'condition2vs1_ttest')
    if any(contains(fieldnames(PE_parameters),'condition'))
        
        condition = PE_parameters.condition;
        if ~isequal(length(condition), size(dataArray,3))
            error('the condition vector needs to be same length as the size(dataArray,3)')
        end
        if ~(  sum(condition==1)==1  && sum(condition==2)==1  )
            error('the conditions vector should include only one "1" and one "2" to compare 2 vs 1')
        end

        nConditions = nnz(condition);
        
        if sum(condition~=0)~=2
            error('i must compare two conditions. rearrange data array')
        end
    else
        error('need to specify the PE_parameters.condition vector for this analysis')
    end
    if any(contains(fieldnames(PE_parameters),'group'))
        group = PE_parameters.group;
        groupLbl = unique(group);
        groupNum = length(groupLbl);
        
        % numerosity of each group
        countTbl = countlabels(group);
        countTbl = sortrows(countTbl,"Label",'ascend'); % enforce ascending order
        if ~isequal(double(string(countTbl.Label)),groupLbl); warning('error with grouping variable'); keyboard; end
    end
end

% correlation
if contains(PE_parameters.test,'correlation')
    corrType = extractAfter(PE_parameters.test,'_');
    if any(contains(fieldnames(PE_parameters),'correlateVar'))
        correlateVar_orig = PE_parameters.correlateVar;
        correlateVar      = PE_parameters.correlateVar;
        if size(dataArray,1)~=size(correlateVar,1)
            error('the var to correlate should have the same num of rows as dataArray')
        end
    else
        error('i need to have a correlation variable. add it to PE_parameters.correlateVar')
    end

    if any(contains(fieldnames(PE_parameters),'group'))
        group = PE_parameters.group;
        groupLbl = unique(group);
        groupNum = length(groupLbl);

        % numerosity of each group
        countTbl = countlabels(group);
        countTbl = sortrows(countTbl,"Label",'ascend'); % enforce ascending order
        if ~isequal(double(string(countTbl.Label)),groupLbl); warning('error with grouping variable'); keyboard; end

    end

    if any(contains(fieldnames(PE_parameters),'condition'))
        condition = PE_parameters.condition;
        nConditions = sum(condition~=0);
    end
end


%% resampling options

rng(PE_parameters.randomSeed)
if strcmp(PE_parameters.stats, 'H0MonteCarlo')

    switch PE_parameters.test

        case 'group2vs1_ttest'
            % shuffle participants with or without replacement
            switch PE_parameters.H0MonteCarlo_replacement
                case false % permutation
                    subjIdxH0 = nan(nSubj,nIterations); % for permutation: shuffle all subjects
                    for itIdx = 1:nIterations
                        subjIdxH0(:,itIdx) = randperm(nSubj,nSubj)';
                    end

                case true % bootstrap
                    subjIdxBoot   = randi(nSubj,nSubj,nIterations); % for bootstrap: sample with replacement
                    subjIdxH0 = nan(nSubj,nIterations); % for bootstrap H0: permute the sample with replacement
                    for itIdx = 1:nIterations
                        subjIdxH0(:,itIdx) = subjIdxBoot(randperm(nSubj),itIdx);  % permuted version under H0
                    end
            end
        
            
        case 'condition2vs1_ttest'

            % permute conditions (swap them randomly)
            % identify which conditions are active in this comparison and will be permuted
            condIdx = 1:length(condition);
            condIdx_inactive   = condition==0;
            condIdx_active     = condition~=0;

            condIdxH0        = repmat(condIdx,nSubj,1,nIterations); % initialize with unchanged order
            for itIdx = 1:nIterations
                for subjIdx = 1:size(dataArray,1)

                    % half of the times, swap the two conditions to compare
                    if rand < 0.5
                        condIdxH0(subjIdx,condIdx_active,itIdx) = flip(condIdx(condIdx_active));
                    end
                     
                end
            end
            
            
            % if replacement it true
            % bootstrap sample of participants with replacement 
            % within their group (groupNum >= 1)
            if PE_parameters.H0MonteCarlo_replacement

                subjIdxH0 = nan(nSubj,nIterations); % initialize
                for groupIdx = 1:groupNum
                    % earlier code (can be deleted if the newer works well)
                    %subjIdxH0_perGroup = randi(countTbl.Count(groupIdx),countTbl.Count(groupIdx),nIterations);
                    %subjIdxH0_perGroup = subjIdxH0_perGroup + (find(group==double(string(countTbl.Label(groupIdx))),1,'first')-1); % add the idx where this group starts
                    %subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:) = subjIdxH0_perGroup;

                    % newer code ---
                    groupMask = group==groupLbl(groupIdx);
                    nGroupSubj = sum(groupMask);
                    groupIndices = find(groupMask);
                    % bootstrap within each group separately
                    for itIdx = 1:nIterations
                        bootIndices = groupIndices(randi(nGroupSubj, [nGroupSubj, 1]));
                        subjIdxH0(groupMask,itIdx) = bootIndices;
                    end
                    % ---
                    
                end
                
                % sanity checks: kept the groups separate (participants are  randomly sampled with replacement ONLY within their group)
                for groupIdx = 1:groupNum
                    tmp = subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:);
                    if ~isequal(groupLbl(groupIdx),unique(group(unique(tmp(:)))))
                        warning('bootstrap mixed participants between groups and it shouldn''t have happened')
                        keyboard
                    end
                end
            end
         
        case {'correlation_Pearson'  'correlation_Spearman'  'correlation_Kendall'}

            % shuffle participants with or without replacement
            switch PE_parameters.H0MonteCarlo_replacement
                case false % permutation

                    subjIdxH0 = nan(nSubj,nIterations); % initialize
                    for groupIdx = 1:groupNum
                        
                        subjIdxH0_perGroup = zeros(countTbl.Count(groupIdx),nIterations);
                        for itIdx = 1:nIterations
                            subjIdxH0_perGroup(:,itIdx) = randperm(countTbl.Count(groupIdx),countTbl.Count(groupIdx));
                        end
                        subjIdxH0_perGroup = subjIdxH0_perGroup + (find(group==double(string(countTbl.Label(groupIdx))),1,'first')-1); % add the idx where this group starts
                        subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:) = subjIdxH0_perGroup;
                    end

                    % sanity checks: kept the groups separate (participants are randomly sampled without replacement ONLY within their group)
                    for groupIdx = 1:groupNum
                        tmp = subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:);
                        if ~isequal(groupLbl(groupIdx),unique(group(unique(tmp(:)))))
                            warning('permutation with n groups mixed participants between groups and it shouldn''t have happened')
                            keyboard
                        end
                    end

                case true % bootstrap

                    subjIdxBoot = nan(nSubj,nIterations); % for bootstrap: sample with replacement
                    subjIdxH0   = nan(nSubj,nIterations); % for H0 bootstrap: sample with replacement and shuffle
                    for groupIdx = 1:groupNum
                        
                        subjIdxBoot_perGroup = randi(countTbl.Count(groupIdx),countTbl.Count(groupIdx),nIterations);
                        subjIdxBoot_perGroup = subjIdxBoot_perGroup + (find(group==double(string(countTbl.Label(groupIdx))),1,'first')-1); % add the idx where this group starts
                        subjIdxBoot(group==double(string(countTbl.Label(groupIdx))),:) = subjIdxBoot_perGroup;

                        subjIdxH0_perGroup = zeros(countTbl.Count(groupIdx),nIterations);
                        for itIdx = 1:nIterations
                            subjIdxH0_perGroup(:,itIdx) = subjIdxBoot_perGroup(randperm(countTbl.Count(groupIdx),countTbl.Count(groupIdx)),itIdx);
                        end
                        subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:) = subjIdxH0_perGroup;
                    end

                    % sanity checks: kept the groups separate (participants are randomly sampled with replacement ONLY within their group)
                    for groupIdx = 1:groupNum
                        tmpBoot = subjIdxBoot(group==double(string(countTbl.Label(groupIdx))),:);
                        if ~isequal(groupLbl(groupIdx),unique(group(unique(tmpBoot(:)))))
                            warning('bootstrap mixed participants between groups and it shouldn''t have happened')
                            keyboard
                        end
                        tmpH0 = subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:);
                        if ~isequal(groupLbl(groupIdx),unique(group(unique(tmpH0(:)))))
                            warning('H0 bootstrap mixed participants between groups and it shouldn''t have happened')
                            keyboard
                        end
                    end

                    % sanity check: for each iteration, H0 and H0 boot have the same participants
                    for groupIdx = 1:groupNum
                        tmpBoot = subjIdxBoot(group==double(string(countTbl.Label(groupIdx))),:);
                        tmpH0 = subjIdxH0(group==double(string(countTbl.Label(groupIdx))),:);

                        for itIdx = 1:nIterations
                            if unique(tmpBoot(:,itIdx)) ~= unique(tmpH0(:,itIdx))
                                error('boot and H0boot have sampled different participants')
                            end
                        end

                    end
            end

    end
end

%%  implementation

% TO DO: write my own ttest2 function for speed
% TO DO: write my own ttest function for speed
% TO DO: double-check formulas for Cohen's d (see BHD study)
% TO DO: make own correlation function using linear algebra for speed


switch PE_parameters.test

    case 'group2vs1_ttest'

        % initialize vars
        clusterMetrics  = repmat(struct('id', [], 'size', [], 'mass', []), 1,nConditions,nIterations);
        clusterMatrix = cell(1,nConditions);
        clustIDList = cell(1,nConditions);

        for itIdx = 1:nIterations
            % shuffle group labels (1st dimension) at each iteration
            if strcmp(PE_parameters.stats,'H0MonteCarlo')
                dataArray(:,:,:) = dataArray_orig(subjIdxH0(:,itIdx),:,:);
            end

            % t-test
            [~,p,~,stats] = ttest2( ...
                dataArray(group==2,:,:), ...
                dataArray(group==1,:,:));

            % Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
%             [p, ~, stats] = ranksum( ...
%                 dataArray(group==2,:,:), ...
%                 dataArray(group==1,:,:))

            tvals = stats.tstat;
            pvals = p;
            % Cohen's d computation for independent samples
            n1 = sum(group==1);
            n2 = sum(group==2);
            s1_squared = var(dataArray(group==1,:,:),0,1);
            s2_squared = var(dataArray(group==2,:,:),0,1);
            pooled_sd = sqrt(((n1-1)*s1_squared+(n2-1)*s2_squared)/(n1+n2-2));
            dCohen = (mean(dataArray(group==2,:,:),1) - mean(dataArray(group==1,:,:),1)) ./ pooled_sd;
        
    
            % cluster identification
            for condIdx = 1:nConditions  % separately per each within-subject variable
                statMatrix = tvals(:,:,condIdx); % choose which metric you want for clusters (eg, tvalues, cohensd)
                pvalMatrix = pvals(:,:,condIdx);
                [clusterMatrix{1,condIdx}, clustIDList{1,condIdx}, clusterMetrics(1,condIdx,itIdx)] = PE_Identification_y1x2z3(statMatrix, pvalMatrix, PE_parameters, DimStruct);
            end

            if strcmp(PE_parameters.stats,'H0MonteCarlo')
                progressBar(itIdx,nIterations)
            end
        end

    case 'condition2vs1_ttest'

        % initialize vars
        clusterMetrics  = repmat(struct('id', [], 'size', [], 'mass', []), groupNum,1,nIterations);
        clusterMatrix = cell(groupNum,1);
        clustIDList = cell(groupNum,1);

        for itIdx = 1:nIterations
            
            if strcmp(PE_parameters.stats,'H0MonteCarlo')

                % shuffle condition info (3rd dimension) at each iteration
                for subjIdx = 1:nSubj
    
                    dataArray(subjIdx,:,:) = dataArray_orig(subjIdx,:,condIdxH0(subjIdx,:,itIdx));
                    
                end

                % if bootstrap option chosen, also bootstrap sample of subjects
                if PE_parameters.H0MonteCarlo_replacement
                    dataArray(:,:,:) = dataArray(subjIdxH0(:,itIdx),:,:);
                end
         
            end

            % this analysis occurs between conditions 2 and 1 separately
            % for as many groups defined (groups can be >= 1), but the
            % threshold are combined. 
            
            tvals = zeros(groupNum, size(dataArray,2) );
            pvals = zeros(groupNum, size(dataArray,2) );
            dCohen = zeros(groupNum, size(dataArray,2) );
            for groupIdx = 1:groupNum % separately per each group (group >= 1)

                % t-test
                [~,p,~,stats] = ttest( ...
                    dataArray(group==groupLbl(groupIdx),:,condition==2), ...
                    dataArray(group==groupLbl(groupIdx),:,condition==1) );

                % Wilcoxon Signed-Rank Test
%                 [p, ~, stats] = signrank( ...
%                     dataArray(group==groupLbl(groupIdx),:,condition==2), ...
%                     dataArray(group==groupLbl(groupIdx),:,condition==1) );

                tvals(groupIdx,:) = stats.tstat;
                pvals(groupIdx,:) = p;
                
                dCohen(groupIdx,:) = mean(dataArray(group==groupLbl(groupIdx),:,condition==2) - dataArray(group==groupLbl(groupIdx),:,condition==1),1)  ./  ...
                    std(dataArray(group==groupLbl(groupIdx),:,condition==2) - dataArray(group==groupLbl(groupIdx),:,condition==1),0,1);
    
                % cluster identification
                statMatrix = tvals(groupIdx,:,:);
                pvalMatrix = pvals(groupIdx,:,:);
                if any(isinf(statMatrix))
                    warning('statMatrix contains infinite values')
                    keyboard
                end
                [clusterMatrix{groupIdx,1}, clustIDList{groupIdx,1}, clusterMetrics(groupIdx,1,itIdx)] = PE_Identification_y1x2z3(statMatrix, pvalMatrix, PE_parameters, DimStruct);
            end
            
            if strcmp(PE_parameters.stats,'H0MonteCarlo')
                progressBar(itIdx,nIterations)
            end

        end

    case {'correlation_Pearson'  'correlation_Spearman'  'correlation_Kendall'}

        % initialize vars
        clusterMetrics  = repmat(struct('id', [], 'size', [], 'mass', []), groupNum, nConditions, nIterations);
        clusterMatrix = cell(groupNum,nConditions);
        clustIDList = cell(groupNum,nConditions);

        for itIdx = 1:nIterations
            % shuffle participants (1st dimension) at each iteration
            if strcmp(PE_parameters.stats,'H0MonteCarlo')
                dataArray(:,:,:) = dataArray_orig(subjIdxH0(:,itIdx),:,:);

                % if using bootstrap, make sure the same participants (not shuffled) are also sampled for the variable to correlate
                if PE_parameters.H0MonteCarlo_replacement
                    correlateVar = correlateVar_orig(subjIdxBoot(:,itIdx));
                    if length(unique(correlateVar))==1
                        fprintf('Skipping iteration %d: all values in correlateVar are identical.\n', itIdx);
                        continue % skip this iteration, as no correlation can be computed
                    end
                       
                end
            end

            % test
            rvals = zeros(groupNum,DimStruct.ny1*DimStruct.nx2*DimStruct.nz3,nConditions);
            pvals = zeros(groupNum,DimStruct.ny1*DimStruct.nx2*DimStruct.nz3,nConditions);
            for groupIdx = 1:groupNum % separately per each group (group >= 1)
                for condIdx = find(condition)' % separately per each condition (condition >= 1)


                    [r, p] = corr(dataArray(group==double(string(countTbl.Label(groupIdx))),:,condIdx),correlateVar(group==double(string(countTbl.Label(groupIdx)))),'type',corrType);
                    rvals(groupIdx,:,condIdx) = r;
                    pvals(groupIdx,:,condIdx) = p;

                    % cluster identification
                    statMatrix = rvals(groupIdx,:,condIdx); % choose which metric you want for clusters (eg, tvalues, cohensd)
                    pvalMatrix = pvals(groupIdx,:,condIdx);
                    if any(isnan(statMatrix))
                        warning('statMatrix contains nan values')
                        keyboard
                    end
                    [clusterMatrix{groupIdx,condIdx}, clustIDList{groupIdx,condIdx}, clusterMetrics(groupIdx,condIdx,itIdx)] = PE_Identification_y1x2z3(statMatrix, pvalMatrix, PE_parameters, DimStruct);

                end
            end
            

            if strcmp(PE_parameters.stats,'H0MonteCarlo')
                progressBar(itIdx,nIterations)
            end
        end
end

%%

% combine statistical results
switch PE_parameters.test
    case {'group2vs1_ttest' 'condition2vs1_ttest'}
        results.tvals = tvals;
        results.pvals = pvals;
        results.dCohen = dCohen;
        results.clusterMatrix = clusterMatrix;
        results.clustIDList = clustIDList;
    case {'correlation_Pearson'  'correlation_Spearman'  'correlation_Kendall'}
        results.rvals = rvals;
        results.pvals = pvals;
        results.clusterMatrix = clusterMatrix;
        results.clustIDList = clustIDList;
end

% compute H0 distribution with its thresholds
if strcmp(PE_parameters.stats,'H0MonteCarlo')

    clusterMeasure_lbl = fieldnames(clusterMetrics);
    clusterMeasure_lbl = clusterMeasure_lbl(~strcmp(clusterMeasure_lbl,'id')); % remove the 'id' column
    clusterMeasure_num = length(clusterMeasure_lbl);

    switch PE_parameters.inference
        case 'Tmax'
            % edit clusterMetrics to keep only the largest per each iteration and condition
            switch PE_parameters.test

                case 'group2vs1_ttest'
                    for condIdx = 1:nConditions
                        for itIdx = 1:nIterations
                            for msIdx = 1:clusterMeasure_num
                 
                                clusterValues = [clusterMetrics(1,condIdx,itIdx).(clusterMeasure_lbl{msIdx})];
                                [maxVal, maxIdx] = max(abs(clusterValues));
                                if maxVal==0; keyboard; end
                                clusterMetrics(1,condIdx,itIdx).(clusterMeasure_lbl{msIdx}) = clusterValues(maxIdx);
                            end
                        end
                    end
    
                    % keep most extreme (by removing the non most extreme) per iteration across conditions
                    for itIdx = 1:nIterations
                        for msIdx = 1:clusterMeasure_num
                            clusterValues = [clusterMetrics(1,:,itIdx).(clusterMeasure_lbl{msIdx})];
                            [~, maxIdx] = max(abs(clusterValues));
                            for condIdx = find(logical(~(maxIdx==(1:length(clusterValues)))))
                                clusterMetrics(1,condIdx,itIdx).(clusterMeasure_lbl{msIdx}) = [];
                            end

                            
                        end
                    end

                case 'condition2vs1_ttest'

                    for groupIdx = 1:groupNum
                        for itIdx = 1:nIterations
                            for msIdx = 1:clusterMeasure_num
                                clusterValues = [clusterMetrics(groupIdx,1,itIdx).(clusterMeasure_lbl{msIdx})];
                                [maxVal, maxIdx] = max(abs(clusterValues));
                                if maxVal==0; keyboard; end
                                clusterMetrics(groupIdx,1,itIdx).(clusterMeasure_lbl{msIdx}) = clusterValues(maxIdx);
                            end
                        end
                    end
                    % keep most extreme (by removing the non most extreme) per iteration across groups
                    for itIdx = 1:nIterations
                        for msIdx = 1:clusterMeasure_num
                            clusterValues = [clusterMetrics(:,1,itIdx).(clusterMeasure_lbl{msIdx})];
                            [~, maxIdx] = max(abs(clusterValues));
                            for groupIdx = find(logical(~(maxIdx==(1:length(clusterValues)))))
                                clusterMetrics(groupIdx,1,itIdx).(clusterMeasure_lbl{msIdx}) = [];
                            end
                            
                        end
                    end


                case { 'correlation_Pearson'  'correlation_Spearman'  'correlatio_Kendall' }
                    for groupIdx = 1:groupNum
                        for condIdx = find(condition)'
                            for itIdx = 1:nIterations
                                for msIdx = 1:clusterMeasure_num
                                    [maxVal, maxIdx] = max(abs(clusterMetrics(groupIdx,condIdx,itIdx).(clusterMeasure_lbl{msIdx})));
                                    if maxVal==0; keyboard; end
                                    clusterMetrics(groupIdx,condIdx,itIdx).(clusterMeasure_lbl{msIdx}) = clusterMetrics(groupIdx,condIdx,itIdx).(clusterMeasure_lbl{msIdx})(maxIdx);
                                end
                            end
                        end
                    end
    
                    % keep most extreme (by removing the non most extreme) per iteration across groups and conditions
                    for itIdx = 1:nIterations
                        for msIdx = 1:clusterMeasure_num
                            clusterValues = [clusterMetrics(:,:,itIdx).(clusterMeasure_lbl{msIdx})];
                            [~, maxIdx] = max(abs(clusterValues));

                            nonExtreme = find(logical(~(maxIdx==(1:length(clusterValues)))));
                            [nonExtreme_group, nonExtreme_cond] = ind2sub([groupNum nConditions],nonExtreme);
                            for neIdx = 1:length(nonExtreme_group)
                                clusterMetrics(nonExtreme_group(neIdx),nonExtreme_cond(neIdx),itIdx).(clusterMeasure_lbl{msIdx}) = [];
                            end
                            
                        end
                    end

            end
            
        case 'Tall'
            % keep clusterMetrics as it is because it already contains all clusters
    end
    
    % compute threshold corresponding with p_crit
    clustThreshold = struct();
    tails = ["oneTail" "twoTails"];
    for msIdx = 1:clusterMeasure_num
        for tIdx = 1:length(tails)
            varLbl = [clusterMeasure_lbl{msIdx} '_' tails{tIdx}];
            switch tIdx
                case 1
                    clustThreshold.(varLbl) = prctile(abs([clusterMetrics.(clusterMeasure_lbl{msIdx})]),100-PE_parameters.cluster_p_crit *100);
                case 2
                    
                    clustThreshold.([varLbl '_lower']) = prctile(([clusterMetrics.(clusterMeasure_lbl{msIdx})]),0  +(PE_parameters.cluster_p_crit/2)*100);
                    clustThreshold.([varLbl '_upper']) = prctile(([clusterMetrics.(clusterMeasure_lbl{msIdx})]),100-(PE_parameters.cluster_p_crit/2)*100);
            end
            
        end
    end


    % view the distribution
    doYouWantTheDist = true;
    if doYouWantTheDist
        figure(4); clf;
        fig = gcf; fig.Units = "normalized"; fig.Position = [0.05 0.2 0.9 0.6];


        for msIdx = 1:clusterMeasure_num
            nexttile;
            tIdx = 1; % one tail
            histogram(abs([clusterMetrics.(clusterMeasure_lbl{msIdx})]), 'FaceColor', [1 0 0],'FaceAlpha',0.5);

            hold on

            line(repmat(clustThreshold.([clusterMeasure_lbl{msIdx} '_' tails{tIdx}]),1,2), get(gca,'YLim'), 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 2)
            tIdx = 2; % two tails
            histogram([clusterMetrics.(clusterMeasure_lbl{msIdx})], 'FaceColor', [0 0 1],'FaceAlpha',0.5);
            line(repmat(clustThreshold.([clusterMeasure_lbl{msIdx} '_' tails{tIdx} '_lower']),1,2), get(gca,'YLim'), 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 2)
            line(repmat(clustThreshold.([clusterMeasure_lbl{msIdx} '_' tails{tIdx} '_upper']),1,2), get(gca,'YLim'), 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 2)

            yyaxis right
            nClusters = length([clusterMetrics.(clusterMeasure_lbl{msIdx})]);
            plot([clusterMetrics.(clusterMeasure_lbl{msIdx})],randn(1,nClusters)+ones(1,nClusters),'x','Color',[0 0 1]);
            set(gca,'YTick',[])

            legend(["abs(measure)" "threshold_{abs}" "measure" "thresholds"],'Location','Best')

            %xLim = clustThreshold.(clusterMeasure_lbl{msIdx}) .* [-3 3];
            %set(gca,'XLim',xLim)
            title(clusterMeasure_lbl{msIdx})
        end
        switch PE_parameters.test
            case 'group2vs1_ttest'
                tltLbl = {
                    ['inference type:' PE_parameters.inference]
                    ['nIterations: ' num2str(nIterations)]
                    ['conditions pooled: ' num2str(nConditions)]
                    };
            case 'condition2vs1_ttest'
                tltLbl = {
                    ['inference type:' PE_parameters.inference]
                    ['nIterations: ' num2str(nIterations)]
                    ['groups pooled: ' num2str(groupNum)]
                    
                    };
            case {'correlation_Pearson'  'correlation_Spearman'  'correlation_Kendall'}
                tltLbl = {
                    ['inference type:' PE_parameters.inference]
                    ['nIterations: ' num2str(nIterations)]
                    ['groups pooled: ' num2str(groupNum)]
                    ['conditions pooled: ' num2str(nConditions)]
                    };
        end
        sgtitle(tltLbl );
        
        

    end
    
else
    clustThreshold = [];
end


end % ends the main function

%% progress bar local function
function progressBar(itIdx,nIterations)
    if itIdx==1; fprintf(['iterations (of ' num2str(nIterations) '): ']); end
    if any(itIdx==round(logspace(log10(1),log10(nIterations),40)))
        fprintf([num2str(itIdx) ' '])
    end
    if itIdx==nIterations; fprintf(' ...iterations completed \n'); end
end


%                 if itIdx == 1
%                     progressBarLength = 100;  % total length of the bar
%                     fprintf('Iterations (of %d):\n', nIterations);
%                 end
%                 
%                 % update only at selected log-spaced iterations
%                 if any(itIdx == round(logspace(log10(1), log10(nIterations), progressBarLength)))
%                     progressFraction = itIdx / nIterations;
%                     filledLength = round(progressBarLength * progressFraction);
%                     barStr = [repmat('#', 1, filledLength), repmat('-', 1, progressBarLength - filledLength)];
%                     % Use '\r' to return to the beginning of the line, updating in place
%                     fprintf('\r[%s] %d/%d (%.1f%%)', barStr, itIdx, nIterations, progressFraction * 100);
%                 end
%                 
%                 if itIdx == nIterations
%                     fprintf('\nIterations done.\n');
%                 end