function [clusterMatrix, clustIDList, clusterMetrics] = ...
    PCE_Identification_y1x2z3( ...
    statMatrix, ...
    pvalMatrix, ...
    PCE_parameters, ...
    DimStruct)

% This function identifies clusters in three dimensions: y1, x2, z3.
% The nature of the dimensions is described in the DimStruct structure.
% Two examples are below.
%   example #1: y1=Frequency, x2=Time, z3=Channel
%   example #2: y1=Frequency4Amplitude, x2=Frequency4Phase, z3=Channel
% 
%
% INPUT: 
%
% statMatrix:       matrix of statistical values (e.g., rho, tval) 
%                   dimensions: 1=chan, 2=freq, 3=time
% 
% pvalMatrix        matrix of pvalues associated with statMatrix 
%                   dimensions: 1=chan, 2=freq, 3=time
%
% neighborMatrix    matrix (chan x chan format)
%                   1=neighbors, 0=no
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
% p_crit            1 digit represeting critical p value (e.g., 0.05) for clustering purposes only
%
% pixelSign         1 or -1 to search, respectively, for clusters with positive or negative statitical values
%
% distThreshold     Euclidean distance inclusive threshold ( <= of this threshold) to determine extent of the cluster proximity area in the 2d grid
%
% figFlag_cluster   1=plot figures, 0=don't
%
% figFlag_proximity 1=plot figure, 0=don't 
%
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
% Conceptually, imagine a 3d space made of "cubes" or "voxels". Two
% significant voxels belong in the same cluster if they share a whole
% "face". In other words, if they have the same coordinates for two
% dimensions and the other dimension is contiguous or neighboring.
%
% Acknowledgments: 
% inspired by Groppe's find_clusters.m and its subfunction follow_clust.m
% https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox 
% 
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sanity checks

% sanity checks for DimStruct structure
% TO DO

% sanity checks for PCE_parameters structure
% TO DO
% if abs(pixelSign)~=1
%     error('clusterSign must be 1 (positive) or -1 (negative)')
% end
% if ~(figFlag_cluster==0 || figFlag_cluster==1 || figFlag_proximity==0 || figFlag_proximity==1)
%     error('figFlag must be 1=draw figure or 0=don''t')
% end

% check that stat and pval matrices have the same dimensionality
% if ~isequal(size(statMatrix), size(pvalMatrix))
%     error('statMatrix and pvalMatrix must have the same size')
% end

% check that neighborMatrix is square
% if contains(designLabel,'channel','IgnoreCase',true)
%     if ~isequal(size(neighborMatrix,1), size(neighborMatrix,2))
%         error('neighborMatrix must be a square matrix')
%     end
% end

% check that statMatrix first dimension is channels
% if contains(designLabel,'channel','IgnoreCase',true)
%     if ~isequal(size(statMatrix,1), size(neighborMatrix,2))
%         error('size(statMatrix,1) must be same as size(neighborMatrix,1)')
%     end
% end


cLim = [-1 1];

% cluster identification method
methodList = ["threshold" "geometric"];
methodIdx = find(strcmp(methodList,PCE_parameters.method));
if ~isempty(methodIdx)
    if PCE_parameters.verbose
        disp(['using: ' methodList{methodIdx}])
    end
else
    error(['only coded methods are: ' char(join(methodList,', '))])
end

% sanity checks on statMatrix
ny1 = length(DimStruct.y1_vec);
nx2 = length(DimStruct.x2_vec);
nz3 = length(DimStruct.z3_chanlocs);
if ~isequal(size(statMatrix,2),[ny1*nx2*nz3])
    warning('the dimensions of statMatrix do not agree with the y1, x2, z3 DimStruct structure');
    %keyboard
end

% get N-dimensional grids corresponding with the points of all dimensions
[y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1,1:nx2,1:nz3);


distThresholdSquared = DimStruct.euclDistThreshold^2;


% sanity check statMatrix
if any(isinf(statMatrix))
    error('statMatrix contains infinite values')
end
if any(isnan(statMatrix))
    error('statMatrix contains nan values')
end
statMatrix_orig = statMatrix;

%% fig proximity

if PCE_parameters.figFlag_proximity==1
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
    distMatrix_val_thresholded(distMatrix_val <= distThreshold) = 1;
    [~, minIdx] = min(distMatrix_val);
    distMatrix_val_thresholded(minIdx) = 2;
    
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

%% preprocessing

switch methodIdx
    case 1  % EFthreshold

    case 2  % curvatureConsensus
        %% 1d -> 3d

        % reshape from 1d to 3d space
        statMatrix_3d = reshape(statMatrix,[ny1 nx2 nz3]);

        %% gaussian smoothing
        % of 2d continuous dimensions
        sigma_y1 = PCE_parameters.sigma_y1;
        sigma_x2 = PCE_parameters.sigma_x2;

        % gaussian along x2
        if nx2 > 1
            domain_x2 = -ceil(3*sigma_x2):ceil(3*sigma_x2);
        else
            domain_x2 = 1;
        end
        g_x2 = exp(-domain_x2.^2/(2*sigma_x2^2));
        g_x2 = g_x2/sum(g_x2); % normalize gaussian to [0, 1]

        % gaussian along y1
        if ny1 > 1
            domain_y1 = -ceil(3*sigma_y1):ceil(3*sigma_y1);
        else
            domain_y1 = 1;
        end
        g_y1 = exp(-domain_y1.^2/(2*sigma_y1^2));
        g_y1 = g_y1/sum(g_y1); % normalize gaussian to [0, 1]

        % gaussian smoothing one dimension at the time
        for z3Idx = 1:nz3
            temp = conv2(squeeze(statMatrix_3d(:,:,z3Idx)), g_x2, 'same');
            result = conv2(temp, g_y1', 'same');
            statMatrix_3d(:,:,z3Idx) = result;
        end

        % of channel dimension
        % it might not be needed for a low density EEG (eg, 32 channels)
        %... TO DO: find a way to do smoothing in the channel dimension

        %% second partial derivatives

        partialDerivatives = false(1,3); % logic to tell which partial derivatives are computed

        % y1 dimension
        if ny1<=3
            if PCE_parameters.verbose
                warning(['skipping y1 dimension (n = ' num2str(ny1) '). the second derivative is not meaningful with fewer than 3 points. so i will not apply the second derivative on this dimension'])
            end
        else
            partialDerivatives(1) = true;
            d2y1 = zeros(size(statMatrix));
            d2y1(2:end-1,:,:) = statMatrix(3:end,:,:) -2*statMatrix(2:end-1,:,:) +statMatrix(1:end-2,:,:);
            d2y1(1,:,:)       = statMatrix(3,    :,:) -2*statMatrix(2,      :,:) +statMatrix(1,      :,:); % first edge (1), same as value at index 2
            d2y1(end,:,:)     = statMatrix(end,  :,:) -2*statMatrix(end-1,  :,:) +statMatrix(end-2,  :,:); % last edge (end), same as value at index end-1
        end

        % x2 dimension
        if nx2<=3
            if PCE_parameters.verbose
                warning(['skipping x2 dimension (n = ' num2str(nx2) '). the second derivative is not meaningful with fewer than 3 points. so i will not apply the second derivative on this dimension'])
            end
        else
            partialDerivatives(2) = true;
            d2x2 = zeros(size(statMatrix_3d));
            % consider the triplet [a b c] and h1 = distance a to b and h2 = distance c to b
            a = statMatrix_3d(:,1:end-2,:);
            b = statMatrix_3d(:,2:end-1,:);
            c = statMatrix_3d(:,3:end,:);
            h = repmat(diff(DimStruct.x2_vec),ny1,1,nz3);
            h1 = abs(h(:,1:end-1,:));
            h2 = abs(h(:,2:end,:));
            if ~isequal(size(a), size(b), size(c), size(h1), size(h2))
                error('error in the continuous second partial derivative')
            end
            d2x2(:,2:end-1,:) = (h2.*c -(h1+h2).*b +h1.*a) ./ (h1.*h2.*(h1+h2));
            d2x2(:,1,      :) = (h2(:,1,:).*c(:,1,:) -(h1(:,1,:)+h2(:,1,:)).*b(:,1,:) +h1(:,1,:).*a(:,1,:)) ./ (h1(:,1,:).*h2(:,1,:).*(h1(:,1,:)+h2(:,1,:)));
            % isequal(d2x2(:,1,      :), d2x2(:,2,      :)) % first and second must be equal, because we imagine a continuation of the curvature to the edge points
            d2x2(:,end,      :) = (h2(:,end,:).*c(:,end,:) -(h1(:,end,:)+h2(:,end,:)).*b(:,end,:) +h1(:,end,:).*a(:,end,:)) ./ (h1(:,end,:).*h2(:,end,:).*(h1(:,end,:)+h2(:,end,:)));
            % isequal(d2x2(:,end,      :), d2x2(:,end-1,      :)) % last and second last must be equal, because we imagine a continuation of the curvature to the edge points
        end

        % z3 dimension
        if nz3<=20
            if PCE_parameters.verbose
                warning(['skipping y1 dimension (n = ' num2str(ny1) '). the second derivative is not meaningful with fewer than 3 points. so i will not apply the second derivative on this dimension'])
            end
        else
            partialDerivatives(3) = true;
            d2z3 = zeros(size(statMatrix_3d));
            for z3Idx = 1:nz3
                % neighbors: idx, num, and values (in 3d)
                neighIdx = find(DimStruct.z3_Dbinned(z3Idx,:)==1);
                neighNum = length(neighIdx);
                neigh = statMatrix_3d(:,:,neighIdx);
                % neighbors: idx and values (in 3d)
                seedIdx = z3Idx;
                seed = repmat(statMatrix_3d(:,:,seedIdx),1,1,1); % the 3rd dimension does not need to be repeated per neighNum because neighNum in the other variables will be summed across
                % angle between each neighbor and the seed: values (in 3d)
                theta = repmat(reshape(DimStruct.z3_D(neighIdx,seedIdx),1,1,neighNum),ny1,nx2,1);
                if ~isequal(size(neigh),size(theta),[size(seed) neighNum])
                    error('error in the channel second partial derivative')
                end
                d2z3(:,:,z3Idx) = sum(neigh./(theta.^2),3) - seed.*sum(1/(theta.^2),3);
            end
        end

        % sign consensus among the curvature of each dimension
        if partialDerivatives(1)
            if partialDerivatives(2)
                if partialDerivatives(3)
                    % active dimensions [1 1 1]
                    curvatureConsensus = ((sign(d2y1)==sign(d2x2)) & (sign(d2y1)==sign(d2z3))).*sign(d2y1);
                else
                    % active dimensions [1 1 0]
                    curvatureConsensus = (sign(d2y1)==sign(d2x2)).*sign(d2y1);
                end
            else
                if partialDerivatives(3)
                    % active dimensions [1 0 1]
                    error('case not encountered yet. it will be coded at a later time')
                else
                    % active dimensions [1 0 0]
                    error('case not encountered yet. it will be coded at a later time')
                end
            end
        else
            if partialDerivatives(2)
                if partialDerivatives(3)
                    % active dimensions [0 1 1]
                    curvatureConsensus = (sign(d2x2)==sign(d2z3)).*sign(d2x2);
                else
                    % active dimensions [0 1 0]
                    error('case not encountered yet. it will be coded at a later time')
                end
            else
                if partialDerivatives(3)
                    % active dimensions [0 0 1]
                    error('case not encountered yet. it will be coded at a later time')
                else
                    % active dimensions [0 0 0]
                    error('no dimension active')
                end
            end
        end

        %% combine sign of curvature and of values
        curvature = zeros(size(curvatureConsensus));

        % positive curvature & negative value = a negative bump (a valley)
        curvature(curvatureConsensus>0 & statMatrix_3d<0) = -1;
        
        % negative curvature & positive value = a positive bump (a hill)
        curvature(curvatureConsensus<0 & statMatrix_3d>0) = 1;

        clear curvatureConsensus % free space

        %% 3d -> 1d

        % reshape from 1d to 3d space
        curvature = reshape(curvature,[1 ny1*nx2*nz3]);

        clear statMatrix_3d % free space
end

%% cluster forming
% find cluster IDs

adjMatrix = PCE_parameters.adjMatrix;

switch methodIdx
    case 1 % threshold
        searchMatrix = (abs(pvalMatrix)<PCE_parameters.threshold) .* sign(statMatrix);
    case 2 % geometric
        searchMatrix = abs(statMatrix_orig) .* curvature ;
end

clusterMatrix = zeros(size(statMatrix)); % initialize. will show cluster membership
%recruitedMatrix = zeros(size(statMatrix)); % initialize. will show if the point has recruited already  // UNUSED, CAN BE DELETED
queueMatrix = false(size(statMatrix)); % initialize queue matrix
clustID = 0; % initialize / clusters id when multiplied by sign and cluster numerosity at the end of the script

% apply ignore mask
searchMatrix = searchMatrix .* PCE_parameters.ignoreMask;

for pnt_idx = find(searchMatrix)    % loop through the searchMatrix points

    % if this cluster is not assigned, assign it to the next available cluster
    if clusterMatrix(pnt_idx)==0

        pnt_sign = sign(searchMatrix(pnt_idx));
        % the next available positive or negative cluster ID
        if pnt_sign > 0
            clustID = max(clusterMatrix(:)) + pnt_sign;
        else
            clustID = min(clusterMatrix(:)) + pnt_sign;
        end
        clusterMatrix(pnt_idx) = clustID;
        firstRecruiter = true;
        askedMatrix = false(size(statMatrix)); % initialize askedMatrix to keep track of which points have been asked to join this specific cluster
    else
        continue % move on to the next point
    end

    % chase down all neighbors that should belong to this cluster but are not assigned yet
    recruiting = true;
    while recruiting
        
        if firstRecruiter
            queueMatrix(pnt_idx) = true;
            askedMatrix(pnt_idx) = true;
        end

        % find neighbors. points need to be...
        % ...in the adjacency matrix
        criteria = adjMatrix(queueMatrix,:);  % some of these points are counted more than once (same of these points share neighbors)
        criteria = sum(criteria,1)>0;  % this line removed doubles
%         if size(criteria,1)~=1
%             error('criteria should be a row vecor')
%         end

        % ...in the searchMatrix and with the same sign as pnt_idx
        criteria = criteria & sign(searchMatrix)==pnt_sign;

        % ...unassigned to a cluster
        criteria = criteria & clusterMatrix==0;

        % ...not already asked to join this cluster 
        % because of the bulk search, some points are asked multiple times to joint the same cluster
        % TO DO :
        % - create a askedMatrix logic similar to queue, which is reset at the beginning of a new cluster 
        % - do something similar to this: 
        criteria = criteria & ~askedMatrix;


        % assign these points to this cluster
        clusterMatrix(criteria) = clustID;

        % assign this point to the askMatrix
        askedMatrix(criteria) = true;

        % remove the current point(s) from the queue (their job is done)
        %queueMatrix(recruiter_idx) = false; % older
        queueMatrix(queueMatrix) = false;

        % assign these points to the queue to search their neighbors
        queueMatrix(criteria) = true;

        if firstRecruiter==true
            firstRecruiter = false;
        end

        % stop recruiting for this clusterID if the queue is empty
        if sum(queueMatrix)==0
            recruiting = false;
        end
        %end
    end
    
end

% update clusterMatrix
clusterMatrix(isnan(clusterMatrix)) = 0;

clustIDList = unique(clusterMatrix(clusterMatrix~=0));
nClust = length(clustIDList);

%%

%% compute support areas

% nsupportRegion = 0; % initialize / clusters id when multiplied by sign and cluster numerosity at the end of the script
% if clustID>0
% 
%     supportRegionMatrix = nan(size(statMatrix)); % initialize. will show supportRegion membership
%     recruitedMatrix = zeros(size(statMatrix)); % initialize. will show if the point has recruited already
%     
% 
%     for pnt_idx = find(sign(statMatrix)==pixelSign)'    % loop through all points
% 
%         % if this supportRegion is not assigned, assign it to the next available supportRegion
%         if isnan(supportRegionMatrix(pnt_idx))
%         %if isnan(supportRegionMatrix(pnt_idx)) && sign(statMatrix(pnt_idx))==pixelSign
%             nsupportRegion = nsupportRegion + 1;
%             supportRegionMatrix(pnt_idx) = nsupportRegion;
%             firstRecruiter = true;
%         else
%             continue
%         end
% 
%         recruiting = true;
%         while recruiting
% 
%             % chase down all neighbors of neighbors that should belong to this supportRegion but are not assigned yet
%             if firstRecruiter
%                 recruiter_idx_List = pnt_idx;
%             end
% 
%             for recruiter_idx = recruiter_idx_List
% 
%                 % find neighbors. points need to be...
%                 % ...with the correct sign (positive or negative)
% %                 if firstRecruiter==true
%                     criterion1 = sign(statMatrix)==pixelSign;
% %                 end
%                 % ...unassigned
%                 criterion2 = isnan(supportRegionMatrix);
%                 % ...euclidean distance less than or equal to input value
%                 [pnt_z3, pnt_y1, pnt_x2] = ind2sub(size(supportRegionMatrix), recruiter_idx);
%                 %criterion3 = sqrt( (y1Matrix-pnt_y1).^2 + (x2Matrix-pnt_x2).^2 ) <= distThreshold;
%                 criterion3 = (y1Matrix-pnt_y1).^2 + (x2Matrix-pnt_x2).^2 <= distThresholdSquared;
%                 % ...angular distance less than or equal to input value
%                 % TO DO when also using channels
%                 criterion4 = ones(size(clusterMatrix));
%                 criteria = criterion1 & criterion2 & criterion3 & criterion4;
% 
%                 % assign these points to this cluster
%                 supportRegionMatrix(criteria) = nsupportRegion;
% 
%                 % append the new recruiter to the recruiter list
%                 recruiter_idx_List = [recruiter_idx_List   find(criteria)'];
%                 if firstRecruiter==true
%                     firstRecruiter = false;
%                 end
% 
%                 % remove the current recruiter from the recruiter list
%                 % this recruiter has recruited
%                 recruitedMatrix(recruiter_idx) = 1;
%                 %recruiter_idx_List = setdiff(recruiter_idx_List,find(recruitedMatrix==1));
%                 isRecruited = ismember(recruiter_idx_List,find(recruitedMatrix==1));
%                 recruiter_idx_List = recruiter_idx_List(~isRecruited);
% 
%                 % stop recruiting if there are no new recruiters
%                 if isempty(recruiter_idx_List)
%                     recruiting = false;
%                 end
%             end
%         end
% 
%     end
% 
%     % update clusterMatrix
%     supportRegionMatrix(isnan(clusterMatrix)) = 0;
%     supportRegionMatrix = supportRegionMatrix*pixelSign;  % if supportRegion are negative, make the registration number negative too
% 
% end

%% visual inspection

% active dimensions [0 1 1]
% one 2d plot, % to do: include cluster metric with the contour line

if PCE_parameters.figFlag_cluster==1
    
    % one 2d plot
    if ny1>1  &  nx2>1
        figure(6); clf
        lyt = tiledlayout('flow');
        z3Idx = randi(nz3,1)

        nexttile(lyt)
        mat2plot = squeeze(statMatrix_orig(:,:,z3Idx));
        contour2plot_neg = logical(squeeze(clusterMatrix(z3Idx,:,:)<0));
        contour2plot_pos = logical(squeeze(clusterMatrix(z3Idx,:,:)>0));
        im = imagesc(DimStruct.x2_vec,DimStruct.y1_vec,mat2plot);

        im.Parent.YAxis.Label.String = DimStruct.y1_lbl;
        im.Parent.YScale = 'log';
        im.Parent.YDir = 'normal';
        im.Parent.YLim = DimStruct.y1_vec([1 end]);
        im.Parent.XAxis.Label.String = DimStruct.x2_lbl;
        im.Parent.XLim = DimStruct.x2_vec([1 end]);
        im.Parent.CLim = cLim;

        colormap("turbo")
        hold on
        contour(DimStruct.x2_vec,DimStruct.y1_vec,contour2plot_neg, 1,'LineColor',[1 0 0],'LineWidth',1, 'LineStyle','-');
        contour(DimStruct.x2_vec,DimStruct.y1_vec,contour2plot_pos, 1,'LineColor',[0 0 1],'LineWidth',1, 'LineStyle','-');
        title([ DimStruct.z3_lbl ': ' DimStruct.z3_chanlocs(z3Idx).labels])

        colorbar

        sgtitle(['num of clusters: ' num2str(nClust)])
    end

    % one 1d plot
    if (ny1>1 && nx2==1)  ||  (ny1==1 && nx2>1)
        figure(7); clf
        lyt = tiledlayout('flow');
        z3Idx = randi(nz3,1);
        
        nexttile(lyt)
        if (ny1>1 && nx2==1)
            horizAxisVals = DimStruct.y1_vec;
            horizAxisLbl = DimStruct.y1_lbl;
            horizAxisLim = DimStruct.y1_vec([1 end]);
        elseif (ny1==1 && nx2>1)
            horizAxisVals = DimStruct.x2_vec;
            horizAxisLbl = DimStruct.x2_lbl;
            horizAxisLim = DimStruct.x2_vec([1 end]);
            
        end
        % line to plot
        tmp = reshape(statMatrix_orig,[ny1 nx2 nz3]);
        line2plot = tmp(:,:,z3Idx);
        pl = plot(horizAxisVals,line2plot,'.-','Color',[0 0 0]);
        pl.Parent.XAxis.Label.String = horizAxisLbl;
        pl.Parent.XLim = horizAxisLim;
        hold on

        tmp = reshape(clusterMatrix,[ny1 nx2 nz3]);
        line2plot_negIdx = tmp(:,:,z3Idx) < 0;
        line2plot_posIdx = tmp(:,:,z3Idx) > 0;
        %line(repmat(pl.Parent.XLim',1,2), PCE_parameters.threshold*[-1 1; -1 1],'LineStyle','--')
        line(pl.Parent.XLim, zeros(1,2),'LineStyle','--','Color',[0 0 0 0.2])

        plot(horizAxisVals(line2plot_negIdx),line2plot(line2plot_negIdx),'o','Color',[0 0 1]);
        plot(horizAxisVals(line2plot_posIdx),line2plot(line2plot_posIdx),'o','Color',[1 0 0]);
        hold on
        title([ DimStruct.z3_lbl ': ' DimStruct.z3_chanlocs(z3Idx).labels])

        sgtitle(['num of clusters: ' num2str(nClust)])
    end

    % TO DO:
    % scalp maps
%     figure(7); clf
%     y1Idx = randi(ny1);  % pick a time point at random
%     x2Idx = randi(nx2);  % pick a time point at random
%     radius = 0.70; % for plotting reasons only (radius is not in the computations)
%     nexttile()
%     surrogateData = randn(1,nz3);
%     % all chans
%     topoplot(surrogateData, chanlocs, 'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0.5, 'emarker', {'o','k',3,1}); hold on;  % plot all channels
%     clusters = nonzeros(unique(clusterMatrix(:,y1Idx,x2Idx)));
%     tmpCol = lines(length(clusters));
%     counter = 0;
%     for clusterIdx = clusters'
%         counter = counter +1;
%         % cluster chans
%         chan2plot = find(clusterMatrix(:,y1Idx,x2Idx)==clusterIdx);
%         topoplot(surrogateData(chan2plot), chanlocs(chan2plot),  'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0, 'emarker', {'.',tmpCol(counter,:) ,15,5}); hold on; % plot neighbors of seed channels
%     end
%     title({['y1-dim point: ' num2str(y1Idx)]  ['x2-dim point: ' num2str(x2Idx)]})
%     colormap(lines)
end




% active dimensions [1 1 1]
% one 2d plot (take one random z3 point)
% one scalp map (take one random y1 and x2 points)

%% compute cluster statistics (e.g., size, mass, robust mass)
% - cluster size in points
% - cluster mass (classic), sum of coefficients
% - cluster normalized mass, extent is normalized to the dimensional space % TO DO: NB: this requires all points to be of the same sign in statMatrix which is not always the case for curvatureConsensus

clusterMetrics = struct(); % initialize
clusterMetrics_size = zeros(1,nClust);  % TO DO: relabel size to extent
clusterMetrics_mass = zeros(1,nClust);

for clIdx = 1:nClust
    idx = clusterMatrix==clustIDList(clIdx); % idx of pixels belonging to this cluster

    % sanity check: all values within this cluster have the same sign in the searchMatrix, not the statMatrix
    if (sign(clustIDList(clIdx))*sum(sign(searchMatrix(idx))))~=sum(idx(:))
        error('this cluster contains a mix of positive and negative values. likely a bug in the code within or outside this function')
    end

%     % find which support region this cluster belongs to
%     searching = true;
%     while searching
%         for srIdx = 1:nsupportRegion
%             sr_idx = supportRegionMatrix==srIdx*pixelSign; % idx of points belonging to this supportRegion
% 
%             if sum(sr_idx(:)+idx(:)==2)==sum(idx(:)==1)
%                 searching = false;
%                 break
%             else
%                 if sum(sr_idx(:)+idx(:)==2)>0
%                     error('the cluster is not fully contained in the support area. very likely a bug')
%                 end
%             end
%         end
%     end

    clusterMetrics_size(1,clIdx) = sum(idx(:));
    clusterMetrics_mass(1,clIdx) = sum(statMatrix_orig(idx)); % can be positive or negative depending on content (ie, interest in pos or neg clusters)

end



clusterMetrics.id       = clustIDList;
clusterMetrics.size     = clusterMetrics_size;
clusterMetrics.mass     = clusterMetrics_mass;

% sanity check: no NaN, no Inf among the clusterMetrics
% TO DO 

if any(isinf([clusterMetrics.size clusterMetrics.mass]))
    warning('infinite values in the clusterMetrics')
    keyboard
end

if any(isnan([clusterMetrics.size clusterMetrics.mass]))
    warning('nan values in the clusterMetrics')
    keyboard
end






