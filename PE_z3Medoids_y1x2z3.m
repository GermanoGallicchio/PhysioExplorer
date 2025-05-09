function [z3_MedIdx, z3_wMedIdx] = ...
    PCE_z3Medoids_y1x2z3(statMat, clustIdxMat,angularDistMat)

% function to identify the medoid and weighted medoid in a discrete space expressing sensor locations on scalp (and their mass)
%
% INPUT: 
%
% statMat       matrix containing all statistical values
% clustIdxMat      matrix identifying the cluster that we want to describe
%               1=member of the cluster, 0=otherwise
%
% OUTPUT
%
% medoidIdx     index identfyying the medoid of the cluster
% 
% Author: Germano Gallicchio
% germano.gallicchio@gmail.com

%% debugging cell

%clustIdxMat   = idx;  % TO DO: try with sparse matrix
%statMat       = eSizMatrix;
%angularDistMat = DimStruct.z3_D;

%% sanity checks and get some initial data

% clustIdxMat and statMat have the same dimensionality
if ~isequal(size(clustIdxMat),size(statMat))
    error('clustIdxMat and statMat need to have the same dimensionality')
end

% angular distance matrix should be a n-n matrix (n = num of channels in the dataset)
if ~isequal(size(angularDistMat,1),size(angularDistMat,2))
    error('angularDistMat must be a square matrix')
end

% third dimension should be z3
if ~isequal(size(statMat,3),size(angularDistMat,1))
    error('third dimension of both statMat and clustIdxMat should be z3')
end

% numerosity of each dimension

ny1 = size(statMat,1);
nx2 = size(statMat,2);
nz3 = size(statMat,3);

% 3d grid for each of the two dimensions
[y1Matrix, x2Matrix, z3Matrix]       = ndgrid(1:ny1,1:nx2,1:nz3);


%% 

% find channels within the cluster
chansInCluster = unique(z3Matrix(clustIdxMat)); 
nChansInCluster = length(chansInCluster); % num of chans in the cluster

% angular distance matrix, for channels in the cluster
% (i.e., spacing of a set of n points in angular space)
angDistSubsetMat = angularDistMat(chansInCluster,chansInCluster);
% angular distance vector, for each channel in the cluster
angDistSubsetVec = sum(angDistSubsetMat,1);

% find the robust mass taken by each channel
clustStat = statMat(clustIdxMat);  % statistical value for each point within the cluster
clustChan = z3Matrix(clustIdxMat);  % channel corresponding to each point within the cluster
chanRegularMass = nan(1,nChansInCluster);
for chanIdx = 1:nChansInCluster
    chanStat = clustStat(clustChan==chansInCluster(chanIdx)); % statistical value for each point within the cluster also corresponding a specific channel
    chanRegularMass(chanIdx) = sum(chanStat);  % regular mass, equivalent to sum(chanStat)
end
chanMass = chanRegularMass;
% sanity check
if ~isequal(abs(sum(sign(chanMass))),nChansInCluster)
    warning('the numerosity of single-sign masses does not correspond with the numerosity of channels describing this cluster. nothing to worry about if using the geometric approach')
end

% compute weight H (for _weighted_ medoid)
H = normalize(-abs(chanMass), 'range', [min(angDistSubsetVec) max(angDistSubsetVec)]);

% medoid
% the channel within the cluster with the smallest
% angular distance from the other channels within the same cluster
[~, MedIdx] = min(angDistSubsetVec);
z3_MedIdx = chansInCluster(MedIdx);
% note (not a bug): if there are only two channels their distance from each
% other is the same so no winner. the weighted medoid is more protected
% from this occurrence as it biases towards the channel where the
% statistical mass is larger.


% weighted medoid 
% the channel within the cluster with the smallest product of angular distance 
% and the reciprocal of the absolute value of the mass ocucpied by that channel.
% similar to medoid but we introduce a bias towards channels in the cluster that have 
% larger statistical values
[~, wMedIdx] = min(angDistSubsetVec .* H);
z3_wMedIdx = chansInCluster(wMedIdx);
