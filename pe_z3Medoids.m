function [z3_MedIdx, z3_wMedIdx] = ...
    PE_z3Medoids_y1x2z3(PE_parameters, statMat, clustIdxMat,angularDistMat)
% Identifies medoids in sensor space for EEG/MEG data analysis
%
% This function identifies the medoid and weighted medoid in a discrete space 
% representing sensor locations on the scalp. The medoid represents the most 
% centrally located sensor within a cluster, while the weighted medoid also 
% considers the statistical values at each point.
%
%
% INPUT:
%   statMat         - 3D matrix (y1 × x2 × z3) containing statistical values
%   clustIdxMat     - 3D logical matrix (y1 × x2 × z3) identifying cluster membership
%                     1 = member of cluster, 0 = otherwise
%   angularDistMat  - 2D matrix (nChannels × nChannels) containing angular 
%                     distances between sensors
%
%
% OUTPUT
%   z3_MedIdx       - Index identifying the medoid sensor
%   z3_wMedIdx      - Index identifying the weighted medoid sensor
%
% Example:
%   [z3_MedIdx, z3_wMedIdx] = PE_z3Medoids_y1x2z3(statMat, clustIdxMat,angularDistMat);
%
% 
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% input validation

% check if we have all input variables
narginchk(4,4);

% variables must not be empty
if isempty(statMat) || isempty(clustIdxMat) || isempty(angularDistMat)
    error('input variables mst not be empty');
end

% statMat must be numeric
if ~isnumeric(statMat)
    error('statMat must be numeric');
end

% clustIdxMat must be numeric or logical
if ~(isnumeric(clustIdxMat) || islogical(clustIdxMat))
    error('clustIdxMat must be numeric or logical');
end

% angularDistMat must be numeric
if ~isnumeric(angularDistMat)
    error('angularDistMat must be numeric');
end

% statMat has no NaN or Inf values
if any(isnan(statMat(:))) || any(isinf(statMat(:)))
    error('statMat contains NaN or Inf values');
end

% angularDistMat has no NaN or Inf values
if any(isnan(angularDistMat(:))) || any(isinf(angularDistMat(:)))
    error('angularDistMat contains NaN or Inf values');
end

% sanity check clustIdxMat contains only 0s and 1s
if ~all(ismember(clustIdxMat(:), [0 1]))
    error('clustIdxMat must contain only 0s and 1s');
end

% clustIdxMat and statMat have the same dimensionality
if ~isequal(size(clustIdxMat),size(statMat))
    error('clustIdxMat and statMat must have the same dimensionality')
end

% angular distance matrix should be a n-n matrix (n = num of channels in the dataset)
if ~isequal(size(angularDistMat,1),size(angularDistMat,2))
    error('angularDistMat must be a square matrix')
end

% third dimension should be z3 / third dimension of statMat must match angularDistMat dimensions
if ~isequal(size(statMat,3),size(angularDistMat,1))
    error('third dimension of both statMat and clustIdxMat should be z3')
end

%% initial data

% numerosity of each dimension
[ny1, nx2, nz3] = size(statMat);

% 3d grid for each of the two dimensions
[y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1, 1:nx2, 1:nz3);

% find channels within the cluster
chansInCluster = unique(z3Matrix(clustIdxMat)); 
nChansInCluster = length(chansInCluster); % num of chans in the cluster

%% shortcuts

% handle empty cluster case
if PE_parameters.verbose
    if nChansInCluster == 0
        warning('No channels found in cluster');
        z3_MedIdx = [];
        z3_wMedIdx = [];
        return;
    end
end

% handle single channel case
if nChansInCluster == 1
    z3_MedIdx = chansInCluster;
    z3_wMedIdx = chansInCluster;
    return;
end

%% implementation

% subset of angular distances for channels in cluster
% (i.e., spacing of a set of n points in angular space)
angDistSubsetMat = angularDistMat(chansInCluster,chansInCluster);
% sum of angular distances, for each channel in the cluster
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
% sanity check / mainly for debugging
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
