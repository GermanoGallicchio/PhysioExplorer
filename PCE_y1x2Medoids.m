function [y1_MedIdx, x2_MedIdx, y1_wMedIdx, x2_wMedIdx ] = ...
    ClusterAnalysis_y1x2Medoids(statMat, clustIdxMat)

% function to identify the medoid and weighted medoid in a 2d continuous
% plane (e.g., time-frequency, time-time, frequency-frequency)
%
% INPUT: 
%
% statMat       matrix containing all statistical values
%               size of matrix is (ny1-nx2) 
%               
% clustIdxMat   matrix identifying the cluster that we want to describe
%               1=member of the cluster, 0=otherwise
%               size of matrix is (ny1-nx2) 
%
%
%
% OUTPUT
%
% y1_MedIdx, x2_MedIdx      indices identfying the medoid in 2d space
% 
% y1_wMedIdx, x2_wMedIdx    indices identfying the weighted medoid in 2d space
%
%
%
% Author: Germano Gallicchio 
% germano.gallicchio@gmail.com


%% debugging cell

%clustIdxMat   = idx;  % TO DO: try with sparse matrix
%statMat       = eSizMatrix;

%% sanity checks and get some initial data

% clustIdxMat and statMat have the same dimensionality
if ~isequal(size(clustIdxMat),size(statMat))
    error('clustIdxMat and statMat need to have the same dimensionality')
end

% numerosity of each dimension
ny1 = size(statMat,1);
nx2 = size(statMat,2);

% 2d grid for each of the two dimensions
[y1Mat,  x2Mat]       = ndgrid(1:ny1, 1:nx2);

% numerosity of the cluster
nCl = nnz(clustIdxMat);


%% 

% statistical values for each point within the cluster
clustStat = statMat(clustIdxMat);  



% y1 and x2-dimension coordinates for each point within the cluster
y1Idx = y1Mat(clustIdxMat);
x2Idx = x2Mat(clustIdxMat);

% Euclidean distance matrix
% (i.e., spacing of a set of n points in Euclidean space)
% TO DO: vectorize the following lines to avoid for loops and achieve faster computation time
euclDistance = nan(nCl,nCl);
for rowIdx = 1:nCl
    for colIdx = 1:nCl
        euclDistance(rowIdx,colIdx) = sqrt(  (y1Idx(rowIdx)-y1Idx(colIdx)).^2  +  (x2Idx(rowIdx)-x2Idx(colIdx)).^2 );
    end
end 


% medoid
% the point within the cluster with the smallest
% Euclidean distance score from the other points within the same cluster
[~, MedIdx] = min(sum(euclDistance,2));
y1_MedIdx = y1Idx(MedIdx);
x2_MedIdx = x2Idx(MedIdx);

% compute weight H (for _weighted_ medoid)
H = normalize(-clustStat, 'range', [min(sum(euclDistance,2)) max(sum(euclDistance,2))]);


% weighted medoid 
% the point within the cluster with the smallest product of Euclidean distance 
% score and the reciprocal of the absolute value of that point.
% similar to medoid but we introduce a bias towards points in the cluster that have 
% larger statistical values
[~, wMedIdx] = min(sum(euclDistance,2) .* H);
y1_wMedIdx = y1Idx(wMedIdx);
x2_wMedIdx = x2Idx(wMedIdx);

