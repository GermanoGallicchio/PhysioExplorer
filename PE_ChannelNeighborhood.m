function [ neighborMatrix, D, Dbinned ] = ...
    PE_ChannelNeighborhood(chanlocs, dist_threshold, angDistFig_flag, topoPlotFig_flag, topoPlotFig_chanLbl, eeglabPath)
% Identifies neighboring EEG channels based on angular distance

%   [neighborMatrix, D, Dbinned] = PE_ChannelNeighborhood(chanlocs, dist_threshold, 
%   angDistFig_flag, topoPlotFig_flag, topoPlotFig_chanLbl, eeglabPath) 
%   determines neighboring channels based on their angular distances on a spherical head model
%
%
%   INPUT:
%
%       chanlocs         - EEGLAB channel location structure (EEG.chanlocs)
%                         Must contain fields: sph_theta, sph_phi, labels
%       dist_threshold   - Angular distance threshold in radians below which
%                         channels are considered neighbors
%                         Recommended values:
%                         - 0.36 for 128-channel systems (Biosemi)
%                         - 0.70 for 32-channel systems (Biosemi)
%       angDistFig_flag  - Flag to plot distance matrix and histogram (1=yes, 0=no)
%       topoPlotFig_flag - Flag to plot topographic maps (1=yes, 0=no)
%       topoPlotFig_chanLbl - Cell array of channel labels for topoplots. seed channel(s) for the topoplots. 
%       eeglabPath      - String containing path to EEGLAB (required for topoplots)
%
%
%
%   OUTPUT:
%
%       neighborMatrix   - N×N binary matrix (N = number of channels)
%                         1 = channels are neighbors, 0 = not neighbors
%       D               - N×N matrix of angular distances in radians
%       Dbinned        - N×N matrix of binned distances where:
%                         0 = same channel
%                         1 = within dist_threshold (neighbors)
%                         2 = within 2*dist_threshold
%                         etc.
%
%   This function is inspired by Groppe's spatial_neighbors.m
%   (https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox)
%   but uses angular instead of linear distance. 
%   This function, like Groppe's, assumes electrodes are placed 
%   on a sphere (i.e., sensor sites have same radius)
%
%   Author: Germano Gallicchio (germano.gallicchio@gmail.com)
%
% TO DO: in future, remove dependency from eeglab's topoplot by creating my own topoplot function

%% sanity checks

% Check number of input arguments
narginchk(2, 6);

% default values for optional inputs
if nargin < 3
    angDistFig_flag = 0;
end
if nargin < 4
    topoPlotFig_flag = 0;
end
if nargin < 5
    topoPlotFig_chanLbl = {};
end
if nargin < 6
    eeglabPath = '';
end


% check topoPlotFig_chanLbl if topoplot is requested
if topoPlotFig_flag
    if ~iscell(topoPlotFig_chanLbl) && ~isstring(topoPlotFig_chanLbl)
        error('PE_ChannelNeighborhood:InvalidChannelLabels', ...
              'topoPlotFig_chanLbl must be a cell array or string array of channel labels.');
    end
    
    % convert to cell array if string array
    if isstring(topoPlotFig_chanLbl)
        topoPlotFig_chanLbl = cellstr(topoPlotFig_chanLbl);
    end
    
    % check channel labels exist in chanlocs
    channel_labels = {chanlocs.labels};
    invalid_channels = setdiff(topoPlotFig_chanLbl, channel_labels);
    if ~isempty(invalid_channels)
        error('PE_ChannelNeighborhood:InvalidChannelLabels', ...
              'The following channel labels were not found in chanlocs: %s', ...
              strjoin(invalid_channels, ', '));
    end
end

% check eeglabPath when topoplot is requested
if topoPlotFig_flag && (~ischar(eeglabPath) && ~isstring(eeglabPath))
    error('PE_ChannelNeighborhood:InvalidEEGLABPath', ...
          'eeglabPath must be a character vector or string.');
end

% check eeglabPath exists when provided
if ~isempty(eeglabPath) && ~exist(eeglabPath, 'dir')
    error('PE_ChannelNeighborhood:InvalidEEGLABPath', ...
          'The specified EEGLAB path does not exist: %s', eeglabPath);
end

% topoplot validation
if topoPlotFig_flag && isempty(eeglabPath)
    error('PE_ChannelNeighborhood:MissingEEGLABPath', ...
          'eeglabPath must be provided when topoPlotFig_flag is true.');
end

%% get info

nChan = length(chanlocs);

%% identify channel neighbors based on eeglab chanlocs

% extract angular coordinates of each electrondes (from eeglab chanlocs)
theta = nan(1,nChan); % initialize
phi   = nan(1,nChan); % initialize
for chanIdx = 1:length(chanlocs)
    theta(chanIdx) = chanlocs(chanIdx).sph_theta;
    phi(chanIdx) = chanlocs(chanIdx).sph_phi;
end

% convert degrees to radians
theta = deg2rad(theta);
phi   = deg2rad(phi);

% compute angular distance between all pairs
nChan = length(chanlocs);
D = zeros(nChan, nChan); % initialize D
for chanAIdx = 1:nChan
    for chanBIdx = 1:nChan

        if chanAIdx==chanBIdx
            angDist = 0;
        else
            angDist = acos( sin(phi(chanAIdx))*sin(phi(chanBIdx))  + cos(phi(chanAIdx))*cos(phi(chanBIdx))*cos(theta(chanAIdx)-theta(chanBIdx)) );
        end
        D(chanAIdx,chanBIdx) = angDist;
    end
end

% visualize angular distance

if angDistFig_flag==1
    figure(2); clf
    f = gcf;
    f.Units = 'normalized'; f.Position = [0 0 0.9 0.7];
    nexttile()
    mt = imagesc(1:nChan, 1:nChan, D);
    mt.Parent.YTick = 1:nChan;
    mt.Parent.YTickLabel = string({chanlocs.labels});
    mt.Parent.YLabel.String = 'channel labels';
    mt.Parent.XTick = 1:nChan;
    mt.Parent.XTickLabel = string({chanlocs.labels});
    mt.Parent.XLabel.String = 'channel labels';
    cb = colorbar;
    cb.Label.String = 'angular distance [rad]';

    nexttile()
    counts = tril(D); counts(counts==0)=NaN; % keep only each inter-channel distance
    hst = histogram(counts(:), 'BinWidth',0.1);
    hst.Parent.YLabel.String = {'num of channel pairs' 'at a certain distance'};
    hst.Parent.XLabel.String = 'angular distance [rad]';

end

% discretize sensor distance and define neighbors
%   0 = same channel 
%   1 = within dist_threshold
%   2 = within dist_threshold*2
%   3 = within dist_threshold*3
%   ...
%   n = within Max(D)
% create edges from 0 to max(D) in steps of dist_threshold, making sure the upper edge is included too
angularEdges = 0:dist_threshold:max(D(:));
if angularEdges(end)~=max(D(:));     angularEdges(end+1) = max(D(:));  end
Dbinned = -999*ones(size(D)); % initialize
dIdx = D==angularEdges(1);
Dbinned(dIdx) = angularEdges(1);
for eIdx = 1:length(angularEdges)-1
    dIdx = D>angularEdges(eIdx)  &  D<=angularEdges(eIdx+1);
    Dbinned(dIdx) = eIdx;
end
%   neighbor = if discretized distance == 1 (a sensor cannot be neighbor with itself)
neighborMatrix = Dbinned==1;
neighborNum = sum(neighborMatrix,1);  % number of neighbors for each channel



% topoplot

% plot channel neighbors
if topoPlotFig_flag==1

    % check EEGLAB path 
    if ~exist('topoplot', 'file')
        if isempty(eeglabPath)
            error('PE_ChannelNeighborhood:EEGLABNotFound', ...
                  'EEGLAB''s topoplot function not found. Valid eeglabPath must be provided.');
        end
        addpath(eeglabPath);
        if ~exist('topoplot', 'file')
            error('PE_ChannelNeighborhood:EEGLABNotFound', ...
                  'topoplot function not found in the provided EEGLAB path.');
        end
    end


    % find idx of seed channels
    chanLbl = topoPlotFig_chanLbl; % seed channels
    chanIdx = nan(1,length(chanLbl));
    for cIdx = 1:length(chanLbl)
        chanIdx(cIdx) = find(strcmp(string({chanlocs.labels}), chanLbl(cIdx)));
    end

    figure(4); clf
    f = gcf;
    f.Units = 'normalized'; f.Position = [0.1 0.05 0.45 0.4];
    
    if length(chanLbl)==8
        nrow = 2;
        ncol = 4;
        tld = tiledlayout(nrow,ncol);
    else
        tld = tiledlayout('flow');
    end
    tld.Padding = 'tight';
    tld.TileSpacing = 'tight';
    
    radius = 0.70; % for plotting reasons only (radius is not in the computations)
    xyLim = [-0.455 0.455];
    plotCounter = 0;
    for cIdx = 1:length(chanLbl)
        plotCounter = plotCounter + 1;
    
        nexttile(tld)
        surrogateData = randn(1,nChan);
        % all chans
        topoplot(surrogateData, chanlocs, 'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0.5, 'emarker', {'o', 0.6*[1 1 1],3,1}); hold on;  % plot all channels
        set(gca,'XLim',xyLim, 'YLim',xyLim)
        % seed chan
        chan2plot = chanIdx(cIdx);
        topoplot(surrogateData(chan2plot ), chanlocs(chan2plot),  'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0, 'emarker', {'.',[0    0.4470    0.7410], 15,5}); hold on; % plot seed channel
        set(gca,'XLim',xyLim, 'YLim',xyLim)
        % neighbor chans
        chan2plot = neighborMatrix(chanIdx(cIdx),:);
        if sum(chan2plot)>=1
            topoplot(surrogateData(chan2plot ), chanlocs(chan2plot),  'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0, 'emarker', {'.',[0.8500    0.3250    0.0980], 15,5}); hold on; % plot neighbors of seed channels
            set(gca,'XLim',xyLim, 'YLim',xyLim)
        end
        %ttl = title({['seed channel: ' chanLbl{cIdx}] [num2str(sum(chan2plot)) ' neighbors' ]});
        ttl = title(['seed: ' chanLbl{cIdx}]);
        ttl.FontWeight = 'normal';
        ttl.VerticalAlignment = 'bottom';
        ttl.FontSize = 10;
        %ttl.Position(2) = 0.55;
        
    end
    %sgtitle('neighbors (in red) of seed channels (in blue)')
    %doYouWantPrint = "no";
    %figureName = [mainFolder 'figures\ERP\chanNeighbors' '_' referenceLbl];
    %if strcmp(doYouWantPrint,"yes")
    %    print('-djpeg','-r1000', figureName);
    %end

end

% min, max, median of neighbors per channel
disp(['on average, channels have a MEDIAN of      : '  num2str(median(neighborNum)) ' neighbors'])
disp(['the channel with fewest neighbors has (MIN): '  num2str(min(neighborNum)) ' neighbors'])
disp(['the channel with most neighbors has (MAX)  : '  num2str(max(neighborNum)) ' neighbors'])