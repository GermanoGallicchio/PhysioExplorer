function [clusterMatrix_corrected, clusterMetrics_corrected] = ...
    PE_Pruning_y1x2z3(clusterMatrix, clusterMetrics, clustThresholdH0, clustIDList, PE_parameters)
% PCE_Pruning_y1x2z3 removes clusters based on a threshold derived under
% the null hypothesis (H0). It removes clusters in 3D data that do not meet
% the significance threshold derived from the null hypothesis (H0) testing
% computed by the previous function. Under threshold clusters are set to
% zero in the matrix and NaN in the metrics).
%
% Syntax:
%   [clusterMatrix_corrected, clusterMetrics_corrected] = ...
%       PE_Pruning_y1x2z3(clusterMatrix, clusterMetrics, clustThresholdH0, clustIDList, PCE_parameters)
%
% --- FROM HERE ---
% INPUT:
%   clusterMatrix         Matrix (e.g., channels x time x frequency) indicating cluster membership 
%                        for each data point. 0 = no cluster, 1 = cluster "1", 2=cluster "2", etc.
%   clusterMetrics        Structure containing metrics for each cluster (e.g., mass, size, etc.).
%                        Each field should be a vector with one value per cluster.
%   clustThresholdH0      Structure with thresholds derived under H0 for each metric, e.g.,
%                        clustThresholdH0.mass_oneTail for one-tailed thresholding.
%   clustIDList           Vector containing the IDs of clusters to consider, corresponding to 
%                        clusterMetrics.
%   PCE_parameters        Structure with parameters controlling the pruning, including:
%                        - clusterMetricChoice: string specifying which metric to use for pruning.
%
% OUTPUT:
%   clusterMatrix_corrected  Corrected cluster matrix with clusters below threshold removed (zeroed).
%   clusterMetrics_corrected Corrected cluster metrics structure with below-threshold removed (set to NaN).
% 
% 
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% input validation

% same number of clusters in clusterMatrix and clusterMeasure
if length(nonzeros(unique(clusterMatrix(:))))~=length(clusterMetrics.id)
    error('clusterMatrix and clusterMass must have the same number of clusters')
end

% nDim = size(clusterMatrix,2); % length of dimension "1" (e.g., time)
% nChan = size(clusterMatrix,1); % length of dimension "channels"

%% implementation

% idx clusters to remove
clusters2remove_idx  = abs(clusterMetrics.(PCE_parameters.clusterMetricChoice)) <= clustThresholdH0.([PCE_parameters.clusterMetricChoice '_oneTail']);

% correct clusterMatrix
clusterMatrix_corrected = clusterMatrix;  % initialize a copy
for clIdx = find(clusters2remove_idx)
    clusterMatrix_corrected(clusterMatrix==clustIDList(clIdx)) = 0;
end

% correct clusterMeasure
clusterMetrics_lbl = string(fieldnames(clusterMetrics));
clusterMetrics_num = length(clusterMetrics_lbl);
clusterMetrics_corrected = clusterMetrics;  % initialize a copy
for csIdx = 1:clusterMetrics_num 
    clusterMetrics_corrected.(clusterMetrics_lbl(csIdx))(clusters2remove_idx) = NaN;
end


disp(['over a total number of ' num2str(length(clusters2remove_idx)) ' clusters: removed ' num2str(sum(clusters2remove_idx)) ', retained ' num2str(sum(~clusters2remove_idx))])