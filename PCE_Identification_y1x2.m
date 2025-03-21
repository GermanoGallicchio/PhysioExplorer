function [clusterMatrix, clusterSize, clusterMass, clusterRobMass] = ...
    PCE_Identification_y1x2(statMatrix, pvalMatrix, DimStruct, p_crit, pixelSign, distThreshold, figFlag_cluster, figFlag_proximity)

% This function identifies clusters in two dimensions: y1, x2
% The nature of the dimensions is described in the DimStruct structure.
% Two examples are below.
%   example #1: y1=Frequency, x2=Time
%   example #2: y1=Frequency4Amplitude, x2=Frequency4Phase
% 
%
% INPUT: 
%
% statMatrix:       matrix of statistical values (e.g., rho, tval) 
%                   ex dimensions: 1=freq, 2=time
% 
% pvalMatrix        matrix of pvalues associated with statMatrix 
%                   ex dimensions: 1=freq, 2=time
%
%
% DimStruct         structure defining the dimensions. Following example #1 above
%                     DimStruct.y1_lbl      = 'Freq';  % label for dimension y1
%                     DimStruct.y1_vec      = freqVec; % vector of values. needed only if contFlat==1
%                     DimStruct.y1_units    = 'Hz';    % label for the units of dimension y1. needed only if contFlat==1
%                     DimStruct.x2_lbl      = 'Time';   % label for dimension x2
%                     DimStruct.x2_vec      = timeVec;  % vector of values. needed only if contFlat==1
%                     DimStruct.x2_units    = 's';      % label for the units of dimension x2. needed only if contFlat==1
%
% p_crit            1 digit represeting critical p value (e.g., 0.05) for clustering purposes only
%
% pixelSign         1 or -1 to search, respectively, for clusters with positive or negative statitical values
%
% figFlag           1=plot figure, 0=don't
%
% OUTPUT:
%
% clusterMatrix 
%
% clusterSize 
% 
% clusterMass 
% 
% clusterRobMass
%
% 
% Lay-term exaplanation: 
% Conceptually, imagine a 2d space made of "squares" or "pixels". Two
% significant pixels belong in the same cluster if they share a whole
% "edge". In other words, if they have the same coordinates for one
% dimension and the other dimension is contiguous.
%
% Acknowledgments: 
% inspired by Groppe's find_clusters.m and its subfunction follow_clust.m
% https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox 
% 
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sanity checks

% check that stat matrix has the expected dimensions
dimVec_obs = size(statMatrix);  % numerosity of each dimension (observed in statMatrix)
dimVec_exp = [length(DimStruct.y1_vec)  length(DimStruct.x2_vec)  ]; % numerosity of each dimension (expected from DimStruct)
if ~isequal(dimVec_exp, dimVec_obs)
    error('dimensionality incongruence between statMatrix and what expected based on DimStruct')
else
    dimVec = dimVec_obs;
    clear dimVec_obs dimVec_exp
end
ny1 = dimVec(1);
nx2 = dimVec(2);


% check that stat and pval matrices have the same dimensionality
if ~isequal(size(statMatrix), size(pvalMatrix))
    error('statMatrix and pvalMatrix must have the same size')
end



if abs(pixelSign)~=1
    error('clusterSign must be 1 (positive) or -1 (negative)')
end

if ~(figFlag_cluster==0 | figFlag_cluster==1 | figFlag_proximity==0 | figFlag_proximity==1)
    error('figFlag must be 1=draw figure or 0=don''t')
end

%% fig proximity

if figFlag_proximity==1
    % continuity structure along the two continuous dimensions
    % compute the Euclidean distance of each point in a square matrix from its central point.

    % create a grid representing the indices of the continuous dimensions, with its center being (0,0)
    distMatrix_halfExtent = 3; % how many indices on each side of the central index?
    [distMatrix_y1, distMatrix_x2] = ndgrid(-distMatrix_halfExtent:distMatrix_halfExtent,-distMatrix_halfExtent:distMatrix_halfExtent);

    % compute the Euclidean distances from the central point (0,0)
    distMatrix_val = sqrt(distMatrix_y1.^2 + distMatrix_x2.^2);

    % cluster search extent
    distMatrix_val_thresholded = double(distMatrix_val <= distThreshold);
    distMatrix_val_thresholded(distMatrix_halfExtent+1,distMatrix_halfExtent+1) = 2;

    % plot
    figure(2); clf
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

    %title(sprintf('cluster search area (threshold ≤ %d)', distThreshold), 'FontWeight', 'normal');
end

%% get data

% initialize clusterMatrix
clusterMatrix = zeros(size(statMatrix));

%%

% apply threshold (based on p < p_crit) and separate pos from neg pixels
threshMatrix = pvalMatrix<p_crit;
threshMatrix = threshMatrix .* sign(statMatrix); % distinguish positive and negative clusters by their sign (+ or -)

% find idx of Above-Threshold Pixels (positive or negative, depending on pixelSign)
ATPixels = find(threshMatrix==pixelSign);

% get 2D grids corresponding with the above-treshold pixels
[y1Matrix, x2Matrix] = ndgrid(1:ny1,1:nx2);


%% implementation

% from here below, "pixels" is used for "above-threshold pixel"
nClust = 0; % initialize / clusters numerosity
pixel_registrationStatus  = false(1,length(ATPixels)); % 0=not belonging to any cluster, 1=member of a cluster
pixel_registrationCluster = nan(1,length(ATPixels));   % number corresponding to the cluster they belong to

for pixel_idx = 1:length(ATPixels)    % loop through above-threshold pixels

    % if this cluster is not registered, register it in the next available cluster
    if pixel_registrationStatus(pixel_idx)==false
        nClust = nClust + 1;
        pixel_registrationStatus(pixel_idx)  = true;
        pixel_registrationCluster(pixel_idx) = nClust;
    end

    % search for neighbors
    toggleSwitch = true; % keep searching until this is true
    while toggleSwitch==true
        recruitCounter = 0;

        for recruiter_idx = find(pixel_registrationCluster==nClust)
            for recruitable_idx = find(pixel_registrationStatus==0)

                % y1 (eg, freq) distance between recruiter and recruitable
                dist_y1 = abs(y1Matrix(ATPixels(recruiter_idx)) - y1Matrix(ATPixels(recruitable_idx)));


                % x2 (eg, time) distance between recruiter and recruitable
                dist_x2 = abs(x2Matrix(ATPixels(recruiter_idx)) - x2Matrix(ATPixels(recruitable_idx)));
                
                
                % decide on this voxel recruitment outcome
                % to join the cluster, at least one dist vars needs to be 0 and the other 1
%                 dist_compound = dist_y1 + dist_x2;
%                 if dist_compound<=1  % the voxel joins this cluster
%                     pixel_registrationCluster(recruitable_idx) = nClust;
%                     pixel_registrationStatus(recruitable_idx) = true;
%                     recruitCounter = recruitCounter + 1;
%                 end
                
                % decide on this voxel recruitment outcome
                % to join the cluster, the criterion below needs to be met
                conditionA = sqrt( dist_y1^2 + dist_x2^2 ) <= distThreshold;  % 2d proximity
                if conditionA
                    pixel_registrationCluster(recruitable_idx) = nClust;
                    pixel_registrationStatus(recruitable_idx) = true;
                    recruitCounter = recruitCounter + 1;
                end
            end
        end
        
        % continue as long as there is a new recruit
        if recruitCounter==0
            toggleSwitch = false;
        end

    end
end

% update clusterMatrix
clusterMatrix(ATPixels) = pixel_registrationCluster;
clusterMatrix = clusterMatrix*pixelSign;  % if clusters are negative, make the registration number negative too

%% compute cluster statistics (e.g., size, mass, robust mass)

clusterSize    = nan(1,nClust); % initialize
clusterMass    = nan(1,nClust); % initialize
clusterRobMass = nan(1,nClust); % initialize
for clIdx = 1:nClust
    idx = clusterMatrix==clIdx*pixelSign; % idx of pixels belonging to this cluster

    clusterSize(clIdx)    = sum(idx(:));
    clusterMass(clIdx)    = mean(statMatrix(idx))   * clusterSize(clIdx); % can be positive or negative depending on content (ie, interest in pos or neg clusters)
    clusterRobMass(clIdx) = median(statMatrix(idx)) * clusterSize(clIdx); % can be positive or negative depending on content (ie, interest in pos or neg clusters)

    % sanity check: all stat values within this cluster have the same sign
    if pixelSign*sum(sign(statMatrix(idx)))~=sum(idx(:))
        error('this cluster contains a mix of positive and negative stats values. likely a bug in the code')
    end

end

%% post-processing sanity check
% plot

if figFlag_cluster==1
    
    % y1 x x2 plot (ie, many 2-d plots)
    figure(6); clf
    lyt = tiledlayout('flow');
    nexttile(lyt)
    mat2plot = squeeze(clusterMatrix(:,:));
    imagesc(1:nx2,1:ny1,mat2plot)
    sgtitle(['num of clusters: ' num2str(nClust)])

end
