function [adjMatrix] = pe_adjacency(pe_cfg, figFlag)
% Creates an adjacency vector for 3D physiological data
%
% generates a sparse adjacency vector describing neighborhood relationships
% across three vectorized dimensions: y1, x2, and z3.
% 
% INPUT: 
%       pe_cfg
%
% OUTPUT:
%       adjMatrix - Sparse N×N adjacency matrix where N = ny1*nx2*nz3
%                  Elements are 1 where points are considered "adjacent"
%                  and might join the same cluster, 0 otherwise
%
% Main idea:
%   The function determines adjacency based on:
%   1. Euclidean distance in the y1-x2 plane
%   2. Angular distance between channels in the spatial dimension
%
%   If two points in the ND space have a Euclidean distance <= than
%   criterion 1 and a angular distance <= than criterion 2 then they are
%   adjacent.
%
%   The thresholds for criteria 1 and 2 are decided by the user in:
%   - pe_cfg.clusterParams.distance_y1x2_euclidean   % units: index 
%   - pe_cfg.clusterParams.distance_z3_angular       % units: radians
%
%   Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% input sanity checks

% clusterParams must exist
fieldNeeded = 'clusterParams';
if ~any(contains(fieldnames(pe_cfg),fieldNeeded))
    error(['pe_cfg needs field: ' fieldNeeded])
end


%% shortcuts 

ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;
Nall = ny1*nx2*nz3;

distance_y1x2_euclidean = pe_cfg.clusterParams.distance_y1x2_euclidean^2; % it's faster because it avoid the sqrt() later
distance_z3_angular = pe_cfg.clusterParams.distance_z3_angular;
z3_distanceMatrix = pe_cfg.clusterParams.z3_distanceMatrix;

%%

% get N-dimensional grids corresponding with the points of all dimensions
[y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1,1:nx2,1:nz3);

%% implementation

% loop through each point and find its adjacent spatial-temporal-spectral points
adjMatrix_logical = false(Nall,Nall);% initialize adjMatrix as logical, to be deleted after converting it to sparse
for pIdx = 1:Nall
    
    [pnt_y1, pnt_x2, pnt_z3] = ind2sub([ny1 nx2 nz3], pIdx);

    % find adjacent points of pIdx
    % ...euclidean distance in y1x2 plane less than or equal to input value
    criterion1 = (y1Matrix-pnt_y1).^2 + (x2Matrix-pnt_x2).^2 <= distance_y1x2_euclidean;
    %find(criterion1)'

    % ...angular distance in channel space less than or equal to input value
    chanNeighIdx = find(z3_distanceMatrix(pnt_z3,:) <= distance_z3_angular); % indices of channels that are neighbors of the channel representing pIdx
    criterion2 = ismember(z3Matrix,chanNeighIdx);

    % combine the two criteria. % this will be one row of the sparse adjMatrix;
    adjMatrix_logical(pIdx,:) = reshape(criterion1 & criterion2,[1 ny1*nx2*nz3]);
    
end

adjMatrix = sparse(adjMatrix_logical);
clear adjMatrix_logical % free memory

%% adjacency matrix figure
if pe_cfg.figFlag
    figure()
    spy(adjMatrix)
    xlabel('all dimensions ny1*nx2*nz3')
    ylabel('all dimensions ny1*nx2*nz3')
    title('adjacency matrix in y1x2z3 space')
    f = gcf;
    f.Units = 'normalized';
    f.Position = [0.6 0.6 0.3 0.3];
end


%% fig proximity

if pe_cfg.figFlag
    
    % continuity structure along the two continuous dimensions
    % compute the Euclidean distance of each point in a square matrix from its central point.

    % create a grid representing the indices of the continuous dimensions, with its center being (0,0)
    distMatrix_halfExtent = 3; % how many indices on each side of the central index?
    if ny1>1
        if nx2>1
            % [y1, x2]
            [distMatrix_y1, distMatrix_x2] = ndgrid(-distMatrix_halfExtent:distMatrix_halfExtent,-distMatrix_halfExtent:distMatrix_halfExtent);
        else
            % [y1, 1]
            [distMatrix_y1, distMatrix_x2] = ndgrid(-distMatrix_halfExtent:distMatrix_halfExtent,0);
        end
    else
        if nx2>1
            % [1, x2]
            [distMatrix_y1, distMatrix_x2] = ndgrid(0,-distMatrix_halfExtent:distMatrix_halfExtent);
        else
            % [1, 1]
            [distMatrix_y1, distMatrix_x2] = 0;
        end
    end
    
    % compute the Euclidean distances from the central point (0,0)
    distMatrix_val = sqrt(distMatrix_y1.^2 + distMatrix_x2.^2);

    % cluster search extent
    distMatrix_val_thresholded = zeros(size(distMatrix_val));
    distMatrix_val_thresholded(distMatrix_val <= distance_y1x2_euclidean) = 1;
    [~, minIdx] = min(distMatrix_val(:));
    distMatrix_val_thresholded(minIdx) = 2;
    
    % plot
    figure(); clf
    f = gcf; f.Units = 'normalized'; f.Position = [0.6    0.05   0.24    0.4];
    imagesc(distMatrix_val_thresholded)

    % x axis
    xlabel('dimension 1');
    set(gca, 'XTick', 0.5:1:size(distMatrix_val, 2)+0.5)  % change ticks for the grid
    set(gca, 'XTickLabel', [])
    % y axis
    ylabel('dimension 2');
    set(gca, 'YTick', 0.5:1:size(distMatrix_val, 1)+0.5)  % change ticks for the grid
    set(gca, 'YTickLabel', [])

    % add gridlines that align with cells
    grid on; % Enable grid lines
    set(gca, 'GridColor', [0 0 0], 'GridAlpha', 1, 'LineWidth', 0.5); % Style the grid

    % graphic elements
    axis equal; % Equal scaling for x and y
    axis tight; % Fit axes tightly around the data

    % colormap
    colormap([0.8 0.8 0.8;     0.8500    0.3250    0.0980;      0    0.4470    0.7410]);

    % display values in the grid
    for rowIdx = 1:size(distMatrix_val,1)
        for colIdx = 1:size(distMatrix_val,2)

            % check if the distance is an integer to determine whether using sqrt symbol or not
            switch mod(distMatrix_val(rowIdx,colIdx), 1)
                case 0
                    textStr = sprintf('%d', distMatrix_val(rowIdx, colIdx)); % display as an integer
                otherwise
                    textStr = sprintf('√%d', distMatrix_y1(rowIdx, colIdx)^2 + distMatrix_x2(rowIdx, colIdx)^2); % display with sqrt symbol
            end

            % decide which color to use for the text
            switch distMatrix_val_thresholded(rowIdx,colIdx)
                case 0
                    textCol = [0 0 0];
                case 1
                    textCol = [1 1 1];
            end

            % overlay the formatted text at each grid cell
            text(colIdx, rowIdx, textStr, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', textCol, 'FontSize', 10, 'FontWeight', 'normal');
        end
    end

    sgtitle(sprintf('cluster search area (threshold ≤ %d)', distance_y1x2_euclidean), 'FontWeight', 'normal');

    f2 = gcf;
    f2.Units = 'normalized';
    f2.Position = [0.5 0.05 0.4 0.4];
end
%% 
% other figure
