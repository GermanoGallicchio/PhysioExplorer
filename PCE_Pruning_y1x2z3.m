function [clusterMatrix_corrected, clusterMetrics_corrected] = ...
    PCE_Pruning_y1x2z3(clusterMatrix, clusterMetrics, clustThresholdH0, clustIDList, PCE_parameters)

% ClusterAnalysis_Correction_y1x2z3 removes clusters in the 3 dimensions
% based on the H0 derived threshold computed by the previous function
%
% INPUT: 
%
% clusterMatrix     Matrix (chan x time format) 
%                   Computed by the previous function
%                   It describes cluster membership of each pixel
%                   0=no cluster, 1=cluster "1", 2=cluster "2", etc.
%
% clusterMeasure    vector (1 x nClust)
%                   Computed by the previous function
%                   It describes, for each cluster identified, its mass
%
% clusterThreshold  A number (1 x 1)
%                   Computed by no earlier function, so that
%                   it is used to remove clusters in clusterMatrix that are 
%                   below-threshold
%
%
% OUTPUT:
% clusterMatrix_corrected   Same as before but with some pixels zeroed
%                           (those within a below-threshold cluster)       
%
% clusterMass_corrected     Same as before but with some cluster(s) NaNed
%                           (those clusters below the threshold)
% 
% 
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% debugging cell

%clusterMatrix
%clusterMass
% clusterThreshold = clustMassThreshold1;

%% sanity checks

% same number of clusters in clusterMatrix and clusterMeasure
if length(nonzeros(unique(clusterMatrix(:))))~=length(clusterMetrics.id)
    error('clusterMatrix and clusterMass must have the same number of clusters')
end

% nDim = size(clusterMatrix,2); % length of dimension "1" (e.g., time)
% nChan = size(clusterMatrix,1); % length of dimension "channels"

%% implementation

% idx clusters to remove
clusters2remove_idx  = abs(clusterMetrics.(PCE_parameters.clusterMetricChoice)) <= clustThresholdH0.(PCE_parameters.clusterMetricChoice);

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


disp(['over a total number of ' num2str(length(clusters2remove_idx)) ' cluster, removed: ' num2str(sum(clusters2remove_idx)) ' ' ])