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

if strcmp(PCE_parameters.inference,'Tall')
    if strcmp(PCE_parameters.stats, 'H0MonteCarlo')
        warning('Tall is exploratory. It does not control for FWER whereas Tmax does')
    end
end

% PCE_parameters.stats label
switch PCE_parameters.stats
    case 'H0MonteCarlo'
    case 'observed'
    otherwise
        error('choose either H0MonteCarlo or observed')
end

%% implementation

nSubj = size(dataArray,1);
nConditions = size(dataArray,3);

if contains(PCE_parameters.test,'group')
    if ~isequal(length(group), nSubj)
        error('the group variable needs to be same length as the size(dataArray,1)')
    end
end

% iterations
switch PCE_parameters.stats
    case 'H0MonteCarlo'
        nIterations = PCE_parameters.nIterations;
    case 'observed'
        nIterations = 1; % no iterations
end

% subject indices for resampling options
if strcmp(PCE_parameters.stats, 'H0MonteCarlo')
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
end


%clear clusterMetrics
%clusterMetrics = struct()
clusterMetrics  = repmat(struct('id', [], 'size', [], 'mass', []), nConditions,nIterations);

clusterMatrix = cell(1,nConditions);
clustIDList = cell(1,nConditions);
for itIdx = 1:nIterations

    % shuffle group info (1st dimension) at each iteration
    if strcmp(PCE_parameters.stats,'H0MonteCarlo')
        dataArray(:,:,:) = dataArray(subjIdxH0(:,itIdx),:,:);
    end

    % test
    switch PCE_parameters.test 
        case 'groupBvsA_ttest'
            % TO DO: make own ttest2 function for speed
            [~,p,~,stats] = ttest2( ...
                dataArray(strcmp(group,'B'),:,:), ...
                dataArray(strcmp(group,'A'),:,:));
            tvals = stats.tstat;
            pvals = p;
            dCohen = (mean(dataArray(strcmp(group,'B'),:,:),1) - mean(dataArray(strcmp(group,'A'),:,:,:,:),1)) ./ std(dataArray(:,:,:),0,1);
        case 'condBvsA_ttest'
            % TO DO
        case 'Pearson_correlation'
            % TO DO
    end

    % cluster identification
    for condIdx = 1:nConditions
        statMatrix = tvals(:,:,condIdx);
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
            for condIdx = 1:nConditions
                for itIdx = 1:nIterations
                    for msIdx = 1:clusterMeasure_num
                        [~, maxIdx] = max(abs(clusterMetrics(condIdx,itIdx).(clusterMeasure_lbl{msIdx})));
                        clusterMetrics(condIdx,itIdx).(clusterMeasure_lbl{msIdx}) = clusterMetrics(condIdx,itIdx).(clusterMeasure_lbl{msIdx})(maxIdx);
                    end
                end
            end
        case 'Tall'
            % keep clusterMetrics as it is because it already contains all clusters
    end
    
    % compute threshold corresponding with p_crit
    for msIdx = 1:clusterMeasure_num
        clustThreshold.(clusterMeasure_lbl{msIdx}) = prctile([clusterMetrics.(clusterMeasure_lbl{msIdx})],100-PCE_parameters.cluster_p_crit *100);
    end

    % view the distribution
    doYouWantTheDist = true;
    if doYouWantTheDist
        figure(4); clf;
        fig = gcf; fig.Units = "normalized"; fig.Position = [0.05 0.2 0.9 0.6];


        for msIdx = 1:clusterMeasure_num
            nexttile;
            histogram(abs([clusterMetrics.(clusterMeasure_lbl{msIdx})]), 'FaceColor', [1 0 0],'FaceAlpha',0.5);
            hold on
            histogram([clusterMetrics.(clusterMeasure_lbl{msIdx})], 'FaceColor', [0 0 1],'FaceAlpha',0.5);
            line(repmat(clustThreshold.(clusterMeasure_lbl{msIdx}),1,2), get(gca,'YLim'), 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 2)
            legend(["abs(measure)" "measure" "threshold_{abs}"],'Location','Best')
            xLim = clustThreshold.(clusterMeasure_lbl{msIdx}) .* [-3 3];
            set(gca,'XLim',xLim)
            title(clusterMeasure_lbl{msIdx})
        end
        tltLbl = {
            ['nIterations: ' num2str(nIterations)]
            ['non-EEG conditions pooled: ' num2str(nConditions)]
            ['inference type:' PCE_parameters.inference]
            };
        sgtitle(tltLbl );
        
        

    end
    
else
    clustThreshold = [];
end