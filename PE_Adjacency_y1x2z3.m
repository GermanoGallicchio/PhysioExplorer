function [adjMatrix] = PE_Adjacency_y1x2z3(DimStruct, adjFig)


% This functions creates the adjacency matrix
% 
% INPUT: 
%
% DimStruct         structure defining the dimensions. Following example #1 above
%                     DimStruct.y1_lbl      = 'Freq';  % label for dimension y1
%                     DimStruct.y1_contFlag = 1;       % the variable is on one continuous scale 1=yes, 0=no
%                     DimStruct.y1_vec      = freqVec; % vector of values. needed only if contFlat==1
%                     DimStruct.y1_units    = 'Hz';    % label for the units of dimension y1. needed only if contFlat==1
%                     DimStruct.x2_lbl      = 'Time';   % label for dimension x2
%                     DimStruct.x2_contFlag = 1;        % the variable is on one continuous scale 1=yes, 0=no
%                     DimStruct.x2_vec      = timeVec;  % vector of values. needed only if contFlat==1
%                     DimStruct.x2_units    = 's';      % label for the units of dimension x2. needed only if contFlat==1
%                     DimStruct.z3_lbl      = 'Channel';
%                     DimStruct.z3_contFlag = 0;
%                     DimStruct.z3_chanlocs = EEG.chanlocs; % channel locations structure in the eeglab format
%                     DimStruct.z3_neighborMatrix = neighborMatrix; created by ClusterAnalysis_ChannelNeighborhood.m
%
% distThreshold     Euclidean distance inclusive threshold ( <= of this threshold) to determine extent of the cluster proximity area in the 2d grid
%
% OUTPUT:
%
% 
% sparse matrix N * N telling who is neighbor to whom
%
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sanity checks

% sanity checks for DimStruct structure
% TO DO

%% get data

ny1 = length(DimStruct.y1_vec);
nx2 = length(DimStruct.x2_vec);
nz3 = length(DimStruct.z3_chanlocs);
Nall = ny1*nx2*nz3;

% get N-dimensional grids corresponding with the points of all dimensions
[y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1,1:nx2,1:nz3);

distThresholdSquared = DimStruct.euclDistThreshold^2; % it's faster because it avoid the sqrt() later

%% implementation

% loop through each point and find its adjacent spatial-temporal-spectral points
adjMatrix_logical = false(Nall,Nall);% initialize adjMatrix as logical, to be deleted after converting it to sparse
for pIdx = 1:Nall

    [pnt_y1, pnt_x2, pnt_z3] = ind2sub([ny1 nx2 nz3], pIdx);

    % find adjacent points of pIdx
    % ...euclidean distance in temporal-spectral plane less than or equal to input value
    criterion1 = (y1Matrix-pnt_y1).^2 + (x2Matrix-pnt_x2).^2 <= distThresholdSquared;
    
    % ...angular distance in channel space less than or equal to input value
    chanNeighIdx = find(DimStruct.z3_D(pnt_z3,:) <= DimStruct.angDistThreshold); % indices of channels that are neighbors of the channel representing pIdx
    criterion2 = ismember(z3Matrix,chanNeighIdx);

    % combine the two criteria. % this will be one row of the sparse adjMatrix;
    adjMatrix_logical(pIdx,:) = reshape(criterion1 & criterion2,[1 ny1*nx2*nz3]);
    
end

adjMatrix = sparse(adjMatrix_logical);
clear adjMatrix_logical % free memory

if adjFig
    spy(adjMatrix)
end
