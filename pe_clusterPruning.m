function results = pe_clusterPruning(pe_cfg,results)

% Descr
% prune level-1 results based on clusters
%
% INPUT:
%
%   pe_cfg
%
%
% OUTPUT:
%
% results
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% shortcuts

metrics_lbl = fieldnames(results.resampling.metrics);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id')); % remove the 'id' column
metrics_num = length(metrics_lbl);
%%


% zero clusterMembership for non-significant clusters
for fIdx = 1:metrics_num

    % initialize
    results.clusters.(['clusterMembership_' metrics_lbl{fIdx} '_obs_Corrected']) = results.clusters.clusterMembership_obs;
    results.clusters.(['clustIDList_' metrics_lbl{fIdx} '_obs_Corrected'])       = results.clusters.clustIDList_obs;

    % identify non-significant clusters
    clusters2remove_idx  = abs(results.clusters.metrics_obs.(metrics_lbl{fIdx})) <= results.clusters.inference_maxT.thresholds.(metrics_lbl{fIdx});

    % prune away non significant clusters
    for clIdx = find(clusters2remove_idx)
        % zero features corresponding with non-significant clusters
        idx = results.clusters.clusterMembership_obs==results.clusters.clustIDList_obs(clIdx); % idx of features corresponding with this cluster
        results.clusters.(['clusterMembership_' metrics_lbl{fIdx} '_obs_Corrected'])(idx) = 0;
        results.clusters.(['clustIDList_' metrics_lbl{fIdx} '_obs_Corrected'])(clIdx) = NaN;
    end
end

