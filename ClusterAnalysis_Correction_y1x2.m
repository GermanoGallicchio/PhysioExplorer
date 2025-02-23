function [clusterMatrix_corrected, clusterMeasure_corrected] = ...
    ClusterAnalysis_Correction_y1x2(clusterMatrix, clusterMeasure, clusterThreshold)

% ClusterAnalysis_Correction_y1x2 removes clusters in the 2 dimensions
% based on the H0 derived threshold computed by the previous function
%
% INPUT: 
%
% clusterMatrix     Matrix (y1 x2 format) 
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
if length(nonzeros(unique(clusterMatrix)))~=length(clusterMeasure)
    error('clusterMatrix and clusterMass must have the same number of clusters')
end

% nDim = size(clusterMatrix,2); % length of dimension "1" (e.g., time)
% nChan = size(clusterMatrix,1); % length of dimension "channels"

%%

% idx of clusters to remove
clusters2remove_idx = abs(clusterMeasure) <= clusterThreshold;

% correct clusterMatrix
clusterMatrix_corrected = clusterMatrix;  % initialize a copy
for clIdx = find(clusters2remove_idx)
    clusterMatrix_corrected(clusterMatrix==clIdx) = 0;
end

% correct clusterMeasure
clusterMeasure_corrected = clusterMeasure;  % initialize a copy
clusterMeasure_corrected(clusters2remove_idx) = NaN;

disp(['over a total number of ' num2str(length(clusters2remove_idx)) ' cluster, removed: ' num2str(sum(clusters2remove_idx)) ' ' ])