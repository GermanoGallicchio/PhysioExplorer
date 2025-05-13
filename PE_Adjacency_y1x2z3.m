function [adjMatrix] = PE_Adjacency_y1x2z3(DimStruct, adjFig)
% Creates an adjacency matrix for 3D physiological data
%
%   [adjMatrix] = PE_Adjacency_y1x2z3(DimStruct, adjFig) generates a sparse
%   adjacency matrix describing neighborhood relationships across three
%   dimensions: y1, x2, and z3.
%
% 
% INPUT: 
%
%       DimStruct - Structure containing dimension definitions. y1 and x2
%       are two continuous dimensions such as time and frequency (but also
%       time and time or frequency and frequency); z3 is the channel
%       dimension. any of these dimensions can be empty if working in a
%       subspace. The structure has fields:
%           .y1_lbl        - Label for continuous dimension 1 (e.g., 'Freq')
%           .y1_contFlag   - Flag for continuous scale (1=yes, 0=no)
%           .y1_vec        - Vector of values defining the dimension (e.g., frequency values)
%           .y1_units      - Units (e.g., 'Hz')
%           .x2_lbl        - Label for continuous dimension 2 (e.g., 'Time')
%           .x2_contFlag   - Flag for continuous scale (1=yes, 0=no)
%           .x2_vec        - Vector of values defining the dimension (e.g., time points)
%           .x2_units      - Units (e.g., 's', 'ms')
%           .z3_lbl        - Label for channel dimension
%           .z3_contFlag   - Flag for continuous scale (always 0 for channels)
%           .z3_chanlocs   - Channel locations structure (EEGLAB format)
%           .z3_neighborMatrix - Channel neighborhood matrix. created by PE_ChannelNeighborhood.m
%           .euclDistThreshold - Threshold for Euclidean distance in time-frequency
%           .angDistThreshold  - Threshold for angular distance between channels
%           .z3_D         - Distance matrix for channels
%
%       adjFig - Boolean flag to display spy plot of adjacency matrix
%
%
% OUTPUT:
%       adjMatrix - Sparse N×N adjacency matrix where N = ny1*nx2*nz3
%                  Elements are 1 where points are neighbors, 0 otherwise
%
% Main idea:
%   The function determines adjacency based on:
%   1. Euclidean distance in the temporal-spectral plane
%   2. Angular distance between channels in the spatial dimension
%
%   Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% input sanity checks

% check number of input arguments
narginchk(1, 2);

% default value for adjFig if not provided
if nargin < 2
    adjFig = false;
end

% check DimStruct is a structure
if ~isstruct(DimStruct)
    error('PE_Adjacency_y1x2z3:InvalidInput', ...
        'DimStruct must be a structure.');
end

% required fields in DimStruct
required_fields = {'y1_lbl', 'y1_contFlag', 'y1_vec', 'y1_units', ...
                  'x2_lbl', 'x2_contFlag', 'x2_vec', 'x2_units', ...
                  'z3_lbl', 'z3_contFlag', 'z3_chanlocs', 'z3_neighborMatrix', ...
                  'euclDistThreshold', 'angDistThreshold', 'z3_D'};

% check if there are missing fields
missing_fields = setdiff(required_fields, fieldnames(DimStruct));
if ~isempty(missing_fields)
    error('PE_Adjacency_y1x2z3:MissingFields', ...
          'The following fields are missing in DimStruct: %s', ...
          strjoin(missing_fields, ', '));
end

% check data types and dimensions
if ~ischar(DimStruct.y1_lbl) || ~ischar(DimStruct.x2_lbl) || ~ischar(DimStruct.z3_lbl)
    error('PE_Adjacency_y1x2z3:InvalidLabels', 'Labels must be character arrays.');
end

if ~isnumeric(DimStruct.y1_vec) || ~isvector(DimStruct.y1_vec) || ...
   ~isnumeric(DimStruct.x2_vec) || ~isvector(DimStruct.x2_vec)
    error('PE_Adjacency_y1x2z3:InvalidVectors', ...
        'y1_vec and x2_vec must be numeric vectors.');
end

if ~isnumeric(DimStruct.euclDistThreshold) || ~isscalar(DimStruct.euclDistThreshold) || ...
   DimStruct.euclDistThreshold <= 0
    error('PE_Adjacency_y1x2z3:InvalidThreshold', ...
          'euclDistThreshold must be a positive numeric scalar.');
end

if ~isnumeric(DimStruct.angDistThreshold) || ~isscalar(DimStruct.angDistThreshold) || ...
   DimStruct.angDistThreshold <= 0
    error('PE_Adjacency_y1x2z3:InvalidThreshold', ...
          'angDistThreshold must be a positive numeric scalar.');
end

if ~isnumeric(DimStruct.z3_D) || ~ismatrix(DimStruct.z3_D)
    error('PE_Adjacency_y1x2z3:InvalidDistanceMatrix', ...
          'z3_D must be a numeric matrix.');
end

% check continuous flags
if ~isscalar(DimStruct.y1_contFlag) || ~ismember(DimStruct.y1_contFlag, [0,1]) || ...
   ~isscalar(DimStruct.x2_contFlag) || ~ismember(DimStruct.x2_contFlag, [0,1]) || ...
   ~isscalar(DimStruct.z3_contFlag) || ~ismember(DimStruct.z3_contFlag, [0,1])
    error('PE_Adjacency_y1x2z3:InvalidcontFlags', ...
          'contFlag fields must be scalar values of 0 or 1.');
end

% check adjFig is logical
if ~islogical(adjFig) && ~(isnumeric(adjFig) && ismember(adjFig, [0,1]))
    error('PE_Adjacency_y1x2z3:InvalidAdjFig', ...
          'adjFig must be a logical value or 0/1.');
end


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
