function [clusterMembership, clustIDList, clusterMetrics] = ...
    pe_clusterForming(pe_cfg, statVal, pVal)
% Cluster forming algorithm(s). This function identifies clusters in a 3D
% statistical matrix (y1, x2, z3) using methods (described below):
% - threshold

% This function identifies clusters in three dimensions: y1, x2, z3.
% The nature of the dimensions is described in the pe_cfg.dimensions structure.
% Two examples are below.
%   example #1: y1=Frequency, x2=Time, z3=Channel
%   example #2: y1=Frequency4Amplitude, x2=Frequency4Phase, z3=Channel
% 
%
% INPUT: 
%
%   statMatrix      - Matrix of statistical values (e.g., correlation coefficients, t-values).
%                     Dimensions must correspond to the structure defined in pe_cfg.dimensions.
%
%   pvalMatrix      - Matrix of p-values associated with statMatrix, with matching dimensions.
%
%   PE_parameters  - Structure with parameters controlling the clustering process, including:
%                       .adjMatrix       - Adjacency matrix defining neighborhood relationships
%                       .ignoreMask      - Mask indicating points to exclude from clustering
%                       .figFlag_cluster - Plot clusters (1 = yes, 0 = no)
%                       .figFlag_proximity - Plot cluster proximity area (1 = yes, 0 = no)
%                       .verbose         - Verbosity flag
%                       .method          - Clustering method: 'threshold' or 'geometric'
% fields specific to 'threshold' method
%                       .threshold       - Critical p-value for significance (e.g., 0.05) 
% fields specific to 'geometric' method
%                       .sigma_y1, .sigma_x2 - Smoothing parameters for continuous dimensions
%
%
%   pe_cfg.dimensions       - Structure describing each of the three dimensions, for example:
%                       .y1_lbl, .x2_lbl, .z3_lbl         - Labels for each dimension
%                       .y1_vec, .x2_vec                  - Coordinate vectors (if continuous)
%                       .y1_units, .x2_units              - Units for each dimension
%                       .z3_chanlocs                      - Channel locations structure (EEGLAB format)
%                       .z3_neighborMatrix                - Channel adjacency matrix
%                       .euclDistThreshold                - Euclidean distance threshold for proximity
%
% OUTPUT:
%   clusterMatrix   - Matrix the same size as statMatrix, indicating cluster assignments (0 = not clustered).
%   clustIDList     - List of unique cluster IDs found in clusterMatrix.
%   clusterMetrics  - Struct containing metrics for each cluster:
%                       .id      - Cluster IDs
%                       .size    - Number of points in each cluster
%                       .mass    - Sum of statMatrix values in each cluster
%
% How it works: 
% - treshold 
%
% With this method, the function groups together statistically significant
% points in a three-dimensional space defined by y1, x2, and z3. Groups
% together adjacent points where the p-value is below a specified
% threshold. To determine statistically significance points, it uses p
% values (critical p decided in PE_parameters's appropriate field).
%
% 
% - geometric
%
% ** work in progress **
% WIth this method, the function uses geometric criteria, such as local
% curvature, to form clusters.
%
% General functioning:
%
% Points are considered part of the same cluster if they are adjacent
% according to the adjacency matrix (based on the a priori 3d structure)
% and meet the criteria of the selected method (based on statistical
% results).
%
%
%
% Acknowledgments: Inspired by Groppe's Mass_Univariate_ERP_Toolbox
% specifically find_clusters.m and its subfunction follow_clust.m
% https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox as of early 2025
% 
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)





%% shortcuts

ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;

adjacencyMatrix = pe_cfg.clusterParams.adjacencyMatrix;
pThreshold = pe_cfg.clusterParams.clusterFormingThreshold;
ignoreMask = pe_cfg.R_ignore;

%% sanity checks

% cluster identification
if any(isnan(statVal))
    %warning('statVal contains nan values')
    % keyboard
end

if ~isequal(size(statVal,2),[ny1*nx2*nz3])
    warning('the dimensions of statMatrix do not agree with the y1, x2, z3 pe_cfg.dimensions structure');
    keyboard
end

% sanity check statMatrix
if any(isinf(statVal))
    warning('statMatrix contains infinite values')
    keyboard
end
if any(isnan(statVal))
    %warning('statMatrix contains nan values')
    %keyboard
end
statMatrix_orig = statVal;  % possibly delete

%%
% get N-dimensional grids corresponding with the points of all dimensions
% [y1Matrix, x2Matrix, z3Matrix] = ndgrid(1:ny1,1:nx2,1:nz3);
% the line above is unused

%% cluster forming
% find cluster IDs

searchMatrix = (pVal<pThreshold) .* sign(statVal);

clusterMembership = zeros(size(statVal)); % initialize. will show cluster membership
%recruitedMatrix = zeros(size(statMatrix)); % initialize. will show if the point has recruited already  // UNUSED, CAN BE DELETED
queueMatrix = false(size(statVal)); % initialize queue matrix
clustID = 0; % initialize / clusters id when multiplied by sign and cluster numerosity at the end of the script

% apply ignore mask
searchMatrix = searchMatrix .* ~ignoreMask;


for pnt_idx = find(searchMatrix)    % loop through the searchMatrix points

    % if this cluster is not assigned, assign it to the next available cluster
    if clusterMembership(pnt_idx)==0

        pnt_sign = sign(searchMatrix(pnt_idx));
        % the next available positive or negative cluster ID
        if pnt_sign > 0
            clustID = max(clusterMembership(:)) + pnt_sign;
        else
            clustID = min(clusterMembership(:)) + pnt_sign;
        end
        clusterMembership(pnt_idx) = clustID;
        firstRecruiter = true;
        askedMatrix = false(size(statVal)); % initialize askedMatrix to keep track of which points have been asked to join this specific cluster
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
        criteria = adjacencyMatrix(queueMatrix,:);  % some of these points are counted more than once (some of these points share neighbors)
        criteria = sum(criteria,1)>0;  % this line removed doubles
%         if size(criteria,1)~=1
%             error('criteria should be a row vecor')
%         end

        % ...in the searchMatrix and with the same sign as pnt_idx
        criteria = criteria & sign(searchMatrix)==pnt_sign;

        % ...unassigned to a cluster
        criteria = criteria & clusterMembership==0;

        % ...not already asked to join this cluster 
        % because of the bulk search, some points are asked multiple times to joint the same cluster
        % TO DO :
        % - create a askedMatrix logic similar to queue, which is reset at the beginning of a new cluster id
        % - do something similar to this: 
        % / BOTH DONE NOW
        criteria = criteria & ~askedMatrix;


        % assign these points to this cluster
        clusterMembership(criteria) = clustID;

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
    end
    
end


% update clusterMembership
clusterMembership(isnan(clusterMembership)) = 0;

clustIDList = unique(clusterMembership(clusterMembership~=0));
nClust = length(clustIDList);


%% visual inspection

cLim = [-1 1];

% active dimensions [0 1 1]
% one 2d plot, % to do: include cluster metric with the contour line


if pe_cfg.figFlag  && pe_cfg.nIterations<=2
    
    statVal_matrix = reshape(statVal,[ny1 nx2 nz3]);
    clusterMembership_matrix = reshape(clusterMembership ,[ny1 nx2 nz3]);
    
    % one 2d plot
    if ny1>1  &  nx2>1
        figure(); clf
        lyt = tiledlayout('flow');
        z3Idx = randi(nz3,1);

        nexttile(lyt)
        
        mat2plot = squeeze(statVal_matrix(:,:,z3Idx));

        contour2plot_neg = logical(squeeze(clusterMembership_matrix(:,:,z3Idx)<0));
        contour2plot_pos = logical(squeeze(clusterMembership_matrix(:,:,z3Idx)>0));
        im = imagesc(pe_cfg.dimensions.x2_vec,pe_cfg.dimensions.y1_vec,mat2plot);

        im.Parent.YAxis.Label.String = pe_cfg.dimensions.y1_lbl;
        im.Parent.YScale = 'log';
        im.Parent.YDir = 'normal';
        im.Parent.YLim = pe_cfg.dimensions.y1_vec([1 end]);
        im.Parent.XAxis.Label.String = pe_cfg.dimensions.x2_lbl;
        im.Parent.XLim = pe_cfg.dimensions.x2_vec([1 end]);
        im.Parent.CLim = cLim;

        colormap("turbo")
        hold on
        contour(pe_cfg.dimensions.x2_vec,pe_cfg.dimensions.y1_vec,contour2plot_neg, 1,'LineColor',[1 0 0],'LineWidth',1, 'LineStyle','-');
        contour(pe_cfg.dimensions.x2_vec,pe_cfg.dimensions.y1_vec,contour2plot_pos, 1,'LineColor',[0 0 1],'LineWidth',1, 'LineStyle','-');
        title([ pe_cfg.dimensions.z3_lbl ': ' pe_cfg.dimensions.z3_chanLbl(z3Idx)])

        colorbar

        sgtitle(['num of clusters: ' num2str(nClust)])
    end

    % one 1d plot
    if (ny1>1 && nx2==1)  ||  (ny1==1 && nx2>1)
        figure(); clf
        lyt = tiledlayout('flow');
        z3Idx = randi(nz3,1);
        
        nexttile(lyt)
        if (ny1>1 && nx2==1)
            horizAxisVals = pe_cfg.dimensions.y1_vec;
            horizAxisLbl = pe_cfg.dimensions.y1_lbl;
            horizAxisLim = pe_cfg.dimensions.y1_vec([1 end]);
        elseif (ny1==1 && nx2>1)
            horizAxisVals = pe_cfg.dimensions.x2_vec;
            horizAxisLbl = pe_cfg.dimensions.x2_lbl;
            horizAxisLim = pe_cfg.dimensions.x2_vec([1 end]);
            
        end
        % line to plot
        tmp = reshape(statMatrix_orig,[ny1 nx2 nz3]);
        line2plot = tmp(:,:,z3Idx);
        pl = plot(horizAxisVals,line2plot,'.-','Color',[0 0 0]);
        pl.Parent.XAxis.Label.String = horizAxisLbl;
        pl.Parent.XLim = horizAxisLim;
        hold on

        tmp = reshape(clusterMembership,[ny1 nx2 nz3]);
        line2plot_negIdx = tmp(:,:,z3Idx) < 0;
        line2plot_posIdx = tmp(:,:,z3Idx) > 0;
        %line(repmat(pl.Parent.XLim',1,2), PCE_parameters.threshold*[-1 1; -1 1],'LineStyle','--')
        line(pl.Parent.XLim, zeros(1,2),'LineStyle','--','Color',[0 0 0 0.2])

        plot(horizAxisVals(line2plot_negIdx),line2plot(line2plot_negIdx),'o','Color',[0 0 1]);
        plot(horizAxisVals(line2plot_posIdx),line2plot(line2plot_posIdx),'o','Color',[1 0 0]);
        hold on
        title([ pe_cfg.dimensions.z3_lbl ': ' pe_cfg.dimensions.z3_chanlocs(z3Idx).labels])

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

%% compute cluster statistics (e.g., size, mass)
% - cluster size in points
% - cluster mass (classic), sum of coefficients

clusterMetrics = struct(); % initialize
clusterMetrics_size     = zeros(1,nClust);  % TO DO: relabel size to numerosity
clusterMetrics_mass     = zeros(1,nClust);

for clIdx = 1:nClust
    idx = clusterMembership==clustIDList(clIdx); % idx of pixels belonging to this cluster

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

    % cluster size (aka numerosity)
    clusterMetrics_size(1,clIdx) = sum(idx(:)); % TO DO: change this to numerosity (better term than size)


    % cluster mass
    % clusterMass = clusterSum = clusterMean * clusterNumerosity
    %clusterMetrics_mass(1,clIdx) = sum(statMatrix_orig(idx)); % can be positive or negative depending on content (ie, interest in pos or neg clusters)
    clusterMetrics_mass(1,clIdx) = abs(sum(statMatrix_orig(idx))); % can be only positive (the sign is given by the cluster ID)
   
end



clusterMetrics.id        = clustIDList;
clusterMetrics.size      = clusterMetrics_size;
clusterMetrics.mass      = clusterMetrics_mass;

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
