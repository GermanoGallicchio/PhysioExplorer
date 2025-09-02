function z3_distanceMatrix = pe_z3distance(pe_cfg, eeglabPath)

% Computes angular pairwise distances between sensors. The distance matrix
% is later (pe_adjacency) used to identify which channels are in each channel's
% "neighborhood".
% 
% Optionally, draws a multiplot figure
%   - distance matrix 
%   - histogram
%   - topographic map showing channel neighborhood for user-selected channel(s)
%   - that graph I made for the Alcohol study
%
%   INPUT:
%
%       pe_cfg
%
%   OUTPUT:
%
%       z3_distanceMatrix          - N×N matrix of angular distances in radians (N = num of channels)
%
% Acknowledgments:
%   This function is inspired on Groppe's spatial_neighbors.m
%   (https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox)
%   but uses angular instead of linear distance. 
%   This function, like Groppe's, assumes electrodes are placed 
%   on a sphere (i.e., sensor sites have same radius)
%
%   Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% shortcuts

chanlocs = pe_cfg.dimensions.z3_chanLocs;
nz3 = pe_cfg.dimensions.z3_num;

% angular coordinates (in radians)
theta = deg2rad([chanlocs(:).sph_theta]);
phi   = deg2rad([chanlocs(:).sph_phi]);

distance_z3_angular = pe_cfg.clusterParams.distance_z3_angular;

%% angular distance matrix

% compute angular distance between all pairs
z3_distanceMatrix = zeros(nz3, nz3); % initialize D
for chanAIdx = 1:nz3
    for chanBIdx = 1:nz3
        if chanAIdx==chanBIdx
            angDist = 0;
        else
            angDist = acos( sin(phi(chanAIdx))*sin(phi(chanBIdx))  + cos(phi(chanAIdx))*cos(phi(chanBIdx))*cos(theta(chanAIdx)-theta(chanBIdx)) );
        end
        z3_distanceMatrix(chanAIdx,chanBIdx) = angDist;
    end
end


%% discretized angular distance matrix

% discretize sensor distance and define neighbors
%   0 = same channel 
%   1 = within dist_threshold
%   2 = within dist_threshold*2
%   3 = within dist_threshold*3
%   ...
%   n = within Max(D)
% create edges from 0 to max(D) in steps of dist_threshold, making sure the upper edge is included too
angularEdges = 0:distance_z3_angular:max(z3_distanceMatrix(:));
if angularEdges(end)~=max(z3_distanceMatrix(:));     angularEdges(end+1) = max(z3_distanceMatrix(:));  end
Dbinned = -999*ones(size(z3_distanceMatrix)); % initialize
dIdx = z3_distanceMatrix==angularEdges(1);
Dbinned(dIdx) = angularEdges(1);
for eIdx = 1:length(angularEdges)-1
    dIdx = z3_distanceMatrix>angularEdges(eIdx)  &  z3_distanceMatrix<=angularEdges(eIdx+1);
    Dbinned(dIdx) = eIdx;
end
%   neighbor = if discretized distance == 1 (a sensor cannot be neighbor with itself)
neighborMatrix = Dbinned==1;
neighborNum = sum(neighborMatrix,1);  % number of neighbors for each channel



%% figure

if pe_cfg.figFlag==1

    %% figure setup
    figure(); clf
    f = gcf;
    f.Units = 'normalized'; f.Position = [0.1 0.01 0.8 0.9];

    panelA = uipanel('Position',[0  0.4  1  0.6]);
    panelB = uipanel('Position',[0  0    0.5  0.4]);
    panelC = uipanel('Position',[0.5  0    0.5  0.4]);

    %% figure panel 2:
    % topoplot

    % plot channel neighbors
    if pe_cfg.figFlag

        % choose number of random channels
        % choose 8 channels if there are at least 8 channels
        % otherwise choose ceil of 25%
        if nz3 >= 8
            nRandChans = 8;
        else
            nRandChans = ceil(nz3/4);
        end

        % choose random channel idx
        topoPlotFig_chanIdx = sort(randperm(nz3,nRandChans));
        topoPlotFig_chanLbl = string({chanlocs(topoPlotFig_chanIdx).labels});
        % propose the random selection to the user
        list = string({chanlocs(:).labels});
        [idx,tf] = listdlg('ListString',list,'SelectionMode','multiple','InitialValue',topoPlotFig_chanIdx,'ListSize',[300 200],'PromptString',{'pure visualization purposes:' 'view the neighborhood of selected seed channels' 'display optimized for 8 channels' '(i preselected a few randomly)'});
        if tf==0
            idx=topoPlotFig_chanIdx; warning(['You did not select, so i kept the random selection'])
        else
            % update the selection of channels
            topoPlotFig_chanIdx = idx;
            topoPlotFig_chanLbl = string({chanlocs(idx).labels});
        end
        chanIdx = topoPlotFig_chanIdx;
        chanLbl = topoPlotFig_chanLbl;

        % pick tile layout options
        if length(chanLbl)==8
            nrow = 2;
            ncol = 4;
            tldB = tiledlayout(panelB,nrow,ncol);
        else
            tldB = tiledlayout(panelB,'flow');
        end
        tldB.Padding = 'tight';
        tldB.TileSpacing = 'tight';

        % draw z3 plot
        plotCounter = 0;
        for cIdx = 1:length(chanLbl)
            plotCounter = plotCounter + 1;

            nexttile(tldB)

            % 3d coordinates
            coord_3d = chanlocs; % needs sph_theta and sph_phi
            % z3Values per coordinate
            z3Values = zeros(1,nz3); % initialize to zero (uninvolved channels will stay zero)
            z3Values(chanIdx(cIdx)) = 1; % set seed channel to 1
            z3Values(neighborMatrix(chanIdx(cIdx),:)) = 2; % set neighbor channels to 2
            % plot parameters
            params.projectionType = 'azimuthalEquidistant';   % azimuthalEquidistant | azimuthalConformal | orthographic
            params.drawLines = true;
            params.lineWidth = 1;
            params.lineCol = [0 0 0 1];
            params.chanMarkerSize = 5;
            params.chanMarkerChar = 'o';
            params.chanLbl = false; % true | false
            params.colBar = false; % true | false
            params.colMap = [1 1 1; lines(2)];

            pe_z3Plot(coord_3d,z3Values,params)

            
            %ttl = title({['seed channel: ' chanLbl{cIdx}] [num2str(sum(chan2plot)) ' neighbors' ]});
            ttl = title(['seed: ' chanLbl{cIdx}]);
            ttl.FontWeight = 'normal';
            ttl.VerticalAlignment = 'bottom';
            ttl.FontSize = 10;
            %ttl.Position(2) = 0.55;

        end
        sgtitle('neighbors (in red) of seed channels (in blue)')

    end

    %% figure panel 1: 
    % angular distance matrix

    tldA = tiledlayout(panelA,'flow');
    
    nexttile(tldA)
    mt = imagesc(1:nz3, 1:nz3, z3_distanceMatrix);
    mt.Parent.YTick = 1:nz3;
    mt.Parent.YTickLabel = string({chanlocs.labels});
    mt.Parent.YLabel.String = 'channel labels';
    mt.Parent.XTick = 1:nz3;
    mt.Parent.XTickLabel = string({chanlocs.labels});
    mt.Parent.XLabel.String = 'channel labels';
    mt.Parent.Colormap = flipud(copper);
    cb = colorbar;
    cb.Label.String = 'angular distance [rad]';


    
    nexttile(tldA)
    counts = tril(z3_distanceMatrix); counts(counts==0)=NaN; % keep only each inter-channel distance
    hst = histogram(counts(:), 'BinWidth',0.1);
    hst.Parent.YLabel.String = {'num of channel pairs' 'at a certain distance'};
    hst.Parent.XLabel.String = 'angular distance [rad]';
    xline(distance_z3_angular,'Color',[1 0 0],'LineStyle','--','LineWidth',1)

    %% figure panel 3:
    % which threshold yields the most similar number of neighbors across channels


    tldC = tiledlayout(panelC,'flow');
    dVal = logspace(log10(0.5+1),log10(3+1),400)-1;
    MinNumNeigh = nan(size(dVal));
    MaxNumNeigh = nan(size(dVal));
    StdNumNeigh = nan(size(dVal));
    for nIdx = 1:length(dVal)
        MinNumNeigh(nIdx) = min(sum(z3_distanceMatrix<=dVal(nIdx),2)-1);
        MaxNumNeigh(nIdx) = max(sum(z3_distanceMatrix<=dVal(nIdx),2)-1);
        StdNumNeigh(nIdx) = std(sum(z3_distanceMatrix<=dVal(nIdx),2)-1);
    end

    nexttile(tldC)
    plot(dVal,MinNumNeigh); hold on
    plot(dVal,MaxNumNeigh)

    xTick = logspace(log10(0.5+1),log10(3+1),6)-1;
    set(gca,'XScale','log','XTick',xTick,'XLim',[dVal(1) dVal(end)])
    sgtitle('effect of distance threshold on neighbor count')
    ylabel('number of neighbor per channel')
    xlabel('angular distance threshold [rad]')

    yyaxis right
    plot(dVal,StdNumNeigh)

    xline(distance_z3_angular,'Color',0.5*[1 1 1],'LineStyle','--','LineWidth',0.25)
    legend(["minimum" "maximum" "st. deviation" num2str(distance_z3_angular)])

end




%% neighborhood feedback (some descriptive stats) to user

if pe_cfg.verbose
    % min, max, median of neighbors per channel
    disp(['on average, channels have a MEDIAN of      : '  num2str(median(neighborNum)) ' neighbors'])
    disp(['the channel with fewest neighbors has (MIN): '  num2str(min(neighborNum)) ' neighbors'])
    disp(['the channel with most neighbors has (MAX)  : '  num2str(max(neighborNum)) ' neighbors'])
end