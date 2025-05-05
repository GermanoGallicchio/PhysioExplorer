function [results, clusterMetrics, clustThreshold] = ...
    PCE_Stats_y1x2z3(dataArray, PCE_parameters, DimStruct, group)

% <<Description>>
%
% INPUT: 
%
% statMatrix:       matrix of statistical values (e.g., rho, tval) 
%                   dimensions: 1=chan, 2=time
% 
% pvalMatrix        matrix of pvalues associated with statMatrix 
%                   dimensions: 1=chan, 2=time
%
% neighborMatrix    matrix (chan x chan format)
%                   1=neighbors, 0=no
%
% DimStruct         structure defining the dimensions. Following example #1 above
%                     DimStruct.z1_lbl      = 'Channel';
%                     DimStruct.z1_contFlag = 0;
%                     DimStruct.z1_chanlocs = EEG.chanlocs; % channel locations structure in the eeglab format
%                     DimStruct.z1_neighborMatrix = neighborMatrix; created by ClusterAnalysis_ChannelNeighborhood.m
%                     DimStruct.x2_lbl      = 'Time';   % label for dimension x2
%                     DimStruct.x2_contFlag = 1;        % the variable is on one continuous scale 1=yes, 0=no
%                     DimStruct.x2_vec      = timeVec;  % vector of values. needed only if contFlat==1
%                     DimStruct.x2_units    = 's';      % label for the units of dimension x2. needed only if contFlat==1
%
% p_crit            1 digit represeting critical p value (e.g., 0.05) for clustering purposes only
%
% pixelSign         1 or -1 to search, respectively, for clusters with positive or negative statitical values
%
% figFlag           1=plot figure, 0=don't
%
% OUTPUT:
%
% clusterMatrix 
%
% clusterSize 
% 
% clusterMass 
% 
% clusterRobMass
%
% 
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sanity checks

dataArray_orig  = dataArray;

if strcmp(PCE_parameters.inference,'Tall')
    if strcmp(PCE_parameters.stats, 'H0MonteCarlo')
        warning('Tall is exploratory. It does not control for FWER whereas Tmax does')
    end
end

% PCE_parameters.stats label
switch PCE_parameters.stats
    case {'H0MonteCarlo'  'observed'}
        if PCE_parameters.verbose
            disp(['using: ' PCE_parameters.stats])
        end
    otherwise
        error('choose either H0MonteCarlo or observed')
end



% number of iterations
switch PCE_parameters.stats
    case 'H0MonteCarlo'
        nIterations = PCE_parameters.nIterations;
    case 'observed'
        nIterations = 1; % no iterations
end

%% parameters

nSubj = size(dataArray,1);
nConditions = size(dataArray,3);

if any(contains(fieldnames(PCE_parameters),'group'))
    group = PCE_parameters.group;
    groupLbl = unique(group);
    groupNum = length(groupLbl);
end

if any(contains(fieldnames(PCE_parameters),'condition'))
    condition = PCE_parameters.condition;
    if sum(condition~=0)~=2
        error('i must compare two conditions. rearrange data array')
    end
end

if contains(PCE_parameters.test,'group')
    if ~contains(fieldnames(PCE_parameters),'group')
        error('need to specify the PCE_parameters.group vector for this analysis')
    else
        if ~isequal(length(group), nSubj)
            error('the group variable needs to be same length as the size(dataArray,1)')
        end
    end
end


if contains(PCE_parameters.test,'condition')

    if ~any(contains(fieldnames(PCE_parameters),'condition'))
        error('need to specify the PCE_parameters.condition vector for this analysis')
    else

        if ~isequal(length(condition), size(dataArray,3))
            error('the condition vector needs to be same length as the size(dataArray,3)')
        end
        if ~(  sum(condition==1)==1  && sum(condition==2)==1  )
            error('the conditions vector should include only one "1" and one "2" to compare 2 vs 1')
        end
    end
end


% subject indices for resampling options
rng(PCE_parameters.randomSeed)
if strcmp(PCE_parameters.stats, 'H0MonteCarlo')

    switch PCE_parameters.test
        case 'group2vs1_ttest'
            % shuffle participants with or without replacement
            switch PCE_parameters.H0MonteCarlo_replacement
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
            condIdx = 1:nConditions;
            condIdx_inactive   = condition==0;
            condIdx_active     = condition~=0;

            condIdxH0        = repmat(condIdx,nSubj,1,nIterations); % initialize in the unchanged order
            
            for itIdx = 1:nIterations
                for subjIdx = 1:nSubj

                    % half of the times, swap the two conditions to compare
                    if rand < 0.5
                        condIdxH0(subjIdx,condIdx_active,itIdx) = flip(condIdx(condIdx_active));
                    end
                     
                end
            end
            
            
            % if replacement it true
            % bootstrap sample of participants with replacement 
            % within their group (groupNum >= 1)
            if PCE_parameters.H0MonteCarlo_replacement

                % numerosity of each group
                countTbl = countlabels(group);
                countTbl = sortrows(countTbl,"Label",'ascend'); % enforce ascending order
                if ~isequal(double(countTbl.Label),groupLbl); error('error with grouping variable'); end

                subjIdxH0 = nan(nSubj,nIterations); % initialize
                for groupIdx = 1:groupNum
                    subjIdxH0_perGroup = randi(countTbl.Count(groupIdx),countTbl.Count(groupIdx),nIterations);
                    subjIdxH0_perGroup = subjIdxH0_perGroup + (find(group==double(countTbl.Label(groupIdx)),1,'first')-1); % add the idx where this group starts
                    subjIdxH0(group==double(countTbl.Label(groupIdx)),:) = subjIdxH0_perGroup;
                end
                
                % sanity checks: kept the groups separate (participants are sampled randomly and with replacement only within their group)
                for groupIdx = 1:groupNum
                    tmp = subjIdxH0(group==double(countTbl.Label(groupIdx)),:);
                    if ~isequal(groupIdx,unique(group(unique(tmp(:)))))
                        error('BUG: bootstrap mixed participants between groups and it shouldn''t have happened')
                    end
                end
            end
            
    end
end

%%  implementation

switch PCE_parameters.test
    case 'group2vs1_ttest'

        % initialize vars
        clusterMetrics  = repmat(struct('id', [], 'size', [], 'mass', []), nConditions,nIterations);
        clusterMatrix = cell(1,nConditions);
        clustIDList = cell(1,nConditions);

        for itIdx = 1:nIterations
            % shuffle group info (1st dimension) at each iteration
            if strcmp(PCE_parameters.stats,'H0MonteCarlo')
                dataArray(:,:,:) = dataArray_orig(subjIdxH0(:,itIdx),:,:);
            end

            % test
            % TO DO: write my own ttest2 function for speed
            [~,p,~,stats] = ttest2( ...
                dataArray(group==2,:,:), ...
                dataArray(group==1,:,:));

            tvals = stats.tstat;
            pvals = p;
            dCohen = (mean(dataArray(strcmp(group,'B'),:,:),1) - mean(dataArray(strcmp(group,'A'),:,:),1)) ./ std(dataArray(:,:,:),0,1);
        
    
            % cluster identification
            for condIdx = 1:nConditions  % separately per each within-subject variable
                statMatrix = tvals(:,:,condIdx); % choose which metric you want for clusters (eg, tvalues, cohensd)
                [clusterMatrix{1,condIdx}, clustIDList{1,condIdx}, clusterMetrics(condIdx,itIdx)] = PCE_Identification_y1x2z3(statMatrix, PCE_parameters, DimStruct);
            end

            if strcmp(PCE_parameters.stats,'H0MonteCarlo')
                if itIdx==1; fprintf(['iterations (of ' num2str(nIterations) '): ']); end
                if any(itIdx==round(logspace(log10(1),log10(nIterations),40)))
                    fprintf([num2str(itIdx) ' '])
                end
                if itIdx==nIterations; fprintf('\niterations done \n'); end
            end
        end

    case 'condition2vs1_ttest'

        % initialize vars
        clusterMetrics  = repmat(struct('id', [], 'size', [], 'mass', []), groupNum,nIterations);
        clusterMatrix = cell(groupNum,1);
        clustIDList = cell(groupNum,1);

        for itIdx = 1:nIterations
            
            if strcmp(PCE_parameters.stats,'H0MonteCarlo')
                % shuffle condition info (3rd dimension) at each iteration

                for subjIdx = 1:nSubj
    
                    dataArray(subjIdx,:,:) = dataArray_orig(subjIdx,:,condIdxH0(subjIdx,:,itIdx));
                    
                end

                % if bootstrap option chosen, also bootstrap sample of subjects
                if PCE_parameters.H0MonteCarlo_replacement
                    dataArray(:,:,:) = dataArray(subjIdxH0(:,itIdx),:,:);
                end
         
            end

            % this analysis occurs between conditions 2 and 1 separately
            % for as many groups defined (groups can be >= 1), but the
            % threshold are combined. 
            % note: it might be more appropriate in most cases to
            % keep thresholds separate by running this analysis separately
            % per group
            
            tvals = zeros(groupNum, size(dataArray,2) );
            pvals = zeros(groupNum, size(dataArray,2) );
            dCohen = zeros(groupNum, size(dataArray,2) );
            for groupIdx = 1:groupNum % separately per each group (group >= 1)

                % TO DO: write my own ttest function for speed

                [~,p,~,stats] = ttest( ...
                    dataArray(group==groupLbl(groupIdx),:,condition==2), ...
                    dataArray(group==groupLbl(groupIdx),:,condition==1) );

                tvals(groupIdx,:) = stats.tstat;
                pvals(groupIdx,:) = p;
                % TO DO: check accuracy of within-subj Cohen's d computation
                dCohen(groupIdx,:) = mean(dataArray(group==groupLbl(groupIdx),:,condition==2) - dataArray(group==groupLbl(groupIdx),:,condition==1),1)  ./  ...
                    std(dataArray(group==groupLbl(groupIdx),:,condition==2) - dataArray(group==groupLbl(groupIdx),:,condition==1),0,1);
    
                % cluster identification
                statMatrix = tvals(groupIdx,:,:); % choose which metric you want for clusters (eg, tvalues, cohensd)
                if any(isinf(statMatrix))
                    warning('statMatrix contains infinite values')
                    keyboard
                end
                [clusterMatrix{groupIdx,1}, clustIDList{groupIdx,1}, clusterMetrics(groupIdx,itIdx)] = PCE_Identification_y1x2z3(statMatrix, PCE_parameters, DimStruct);
            end
            
            if strcmp(PCE_parameters.stats,'H0MonteCarlo')
                % TO DO: add this as a function at the bottom of this script
                if itIdx==1; fprintf(['iterations (of ' num2str(nIterations) '): ']); end
                if any(itIdx==round(logspace(log10(1),log10(nIterations),40)))
                    fprintf([num2str(itIdx) ' '])
                end
                if itIdx==nIterations; fprintf('\niterations done \n'); end
            end

        end
    case 'Pearson_correlation'
    % TO DO: make own correlation function using linear algebra for speed
end

%%

% combine statistical results
results.tvals = tvals;
results.pvals = pvals;
results.dCohen = dCohen;
results.clusterMatrix = clusterMatrix;
results.clustIDList = clustIDList;

% compute H0 distribution with its thresholds
if strcmp(PCE_parameters.stats,'H0MonteCarlo')

    clusterMeasure_lbl = fieldnames(clusterMetrics);
    clusterMeasure_lbl = clusterMeasure_lbl(~strcmp(clusterMeasure_lbl,'id')); % remove the 'id' column
    clusterMeasure_num = length(clusterMeasure_lbl);

    switch PCE_parameters.inference
        case 'Tmax'
            % edit clusterMetrics to keep only the largest per each iteration and condition
            switch PCE_parameters.test
                case 'group2vs1_ttest'
                    for condIdx = 1:nConditions
                        for itIdx = 1:nIterations
                            for msIdx = 1:clusterMeasure_num
                                [~, maxIdx] = max(abs(clusterMetrics(condIdx,itIdx).(clusterMeasure_lbl{msIdx})));
                                clusterMetrics(condIdx,itIdx).(clusterMeasure_lbl{msIdx}) = clusterMetrics(condIdx,itIdx).(clusterMeasure_lbl{msIdx})(maxIdx);
                            end
                        end
                    end
                case 'condition2vs1_ttest'
                    for groupIdx = 1:groupNum
                        for itIdx = 1:nIterations
                            for msIdx = 1:clusterMeasure_num
                 
                                [maxVal, maxIdx] = max(abs(clusterMetrics(groupIdx,itIdx).(clusterMeasure_lbl{msIdx})));
                                if maxVal==0; keyboard; end
                                clusterMetrics(groupIdx,itIdx).(clusterMeasure_lbl{msIdx}) = clusterMetrics(groupIdx,itIdx).(clusterMeasure_lbl{msIdx})(maxIdx);
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
                    clustThreshold.(varLbl) = prctile(abs([clusterMetrics.(clusterMeasure_lbl{msIdx})]),100-PCE_parameters.cluster_p_crit *100);
                case 2
                    
                    clustThreshold.([varLbl '_lower']) = prctile(([clusterMetrics.(clusterMeasure_lbl{msIdx})]),100-(PCE_parameters.cluster_p_crit/2)*100);
                    clustThreshold.([varLbl '_upper']) = prctile(([clusterMetrics.(clusterMeasure_lbl{msIdx})]),0  +(PCE_parameters.cluster_p_crit/2)*100);
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
        switch PCE_parameters.test
            case 'group2vs1_ttest'
                tltLbl = {
                    ['nIterations: ' num2str(nIterations)]
                    ['non-EEG conditions pooled: ' num2str(nConditions)]
                    ['inference type:' PCE_parameters.inference]
                    };
            case 'condition2vs1_ttest'
%TO DO: improve the title: there are no conditions pooled but groups
                tltLbl = {
                    ['nIterations: ' num2str(nIterations)]
                    ['non-EEG groups pooled: ' num2str(groupNum)]
                    ['inference type:' PCE_parameters.inference]
                    };
        end
        sgtitle(tltLbl );
        
        

    end
    
else
    clustThreshold = [];
end