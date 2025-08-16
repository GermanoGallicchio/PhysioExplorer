function pe_z3Plot(coord_3d, z3Values, params)

% draws a 2D scalp map from 3D data in sensor-space (e.g., channels)
% 
%   INPUT:
%
%       coord_3d
%           .sph_theta                spherical theta (like in eeglab EEG.chanlocs structure)
%           .sph_phi                  spherical phi (like in eeglab EEG.chanlocs structure)
%
%       params
%           .projectionType
%
%
%   OUTPUT: 
%       Just display the scalp map
%
%
%   Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

nChan = size(coord_3d,2);
fieldNames = fieldnames(params);

if any(strcmp(fieldNames,'projectionType'))
    projectionType = params.projectionType;
else
    projectionType = 'azimuthalEquidistant'; % default
end


if any(strcmp(fieldNames,'drawLines'))
    drawLines = params.drawLines;
else
    drawLines = true; % default
end

if any(strcmp(fieldNames,'lineWidth'))
    lineWidth = params.lineWidth;
else
    lineWidth = 1; % default
end

if any(strcmp(fieldNames,'lineCol'))
    lineCol = params.lineCol;
else
    lineCol = [0 0 0 1]; % default
end


if any(strcmp(fieldNames,'chanMarkerSize'))
    chanMarkerSize  = params.chanMarkerSize;
else
    chanMarkerSize = 5; % default
end


if any(strcmp(fieldNames,'chanMarkerChar'))
    chanMarkerChar = params.chanMarkerChar;
else
    chanMarkerChar = 'o'; % default
end


if any(strcmp(fieldNames,'chanLbl'))
    chanLbl = params.chanLbl;
else
    chanLbl = false; % default
end


if any(strcmp(fieldNames,'colBar'))
colBar = params.colBar;
else
    colBar = false; % default
end


if any(strcmp(fieldNames,'colMap'))
    colMap = params.colMap;
else
    colMap = turbo; % default
end


if any(strcmp(fieldNames,'cLim'))
    cLim = params.cLim;
else
    cLim = prctile(z3Values,[0 100]); % default
end


%% reference coordinates
% for the very top and the very right of the circumference
% very top = nasion
% very right = ca. T8
% include 3d coordinates of nasion and T8 as extra channels
% these coordinates are not associated with data
% they simply work as reference to scale the units and make the circumference pass by

% nasion 
coord_3d(nChan+1).sph_theta = 90;
coord_3d(nChan+1).sph_phi   = 0;
coord_3d(nChan+1).labels    = 'nasionRef';


%% project electrode 3D locations to 2D

coord_2d = pe_z3Projection(coord_3d, projectionType);


% scale 2D coordinates based on the top and right reference
% so that they represent respectively height of 1 and length of 1 
% in a [0, 1] axis space
coord_2d(:,1) = (coord_2d(:,1)-min(coord_2d(:,1)));
coord_2d(:,2) = (coord_2d(:,2)-min(coord_2d(:,2)));
%references length to normalize
yRef = coord_2d(nChan+1,2);
xRef = coord_2d(nChan+1,1)*2;

coord_2d(:,1) = coord_2d(:,1) / (xRef-min(coord_2d(:,1)));
coord_2d(:,2) = coord_2d(:,2) / (yRef-min(coord_2d(:,2)));

%% value normalization for colors

z3Values_norm = (z3Values - cLim(1)) / (cLim(2) - cLim(1));
idx = round(1 + z3Values_norm * (size(colMap, 1) - 1));
idx = max(1, min(size(colMap,1), idx)); % clip to colormap range

%% figure

% channel markers
for chanIdx = 1:nChan
    plot(coord_2d(chanIdx,1), coord_2d(chanIdx,2), ...
        chanMarkerChar, ...
        'MarkerSize', chanMarkerSize, ...
        'MarkerFaceColor',colMap(idx(chanIdx),:), ...
        'MarkerEdgeColor',lineCol(1,1:3));
    hold on
end
%plot(coord_2d(chanIdx+1,1), coord_2d(chanIdx+1,2), 'X', 'MarkerSize', chanMarkerSize*2); % nasion

% channel labels
if chanLbl
    text(coord_2d(1:nChan,1), coord_2d(1:nChan,2),{coord_3d(1:nChan).labels},'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8)
end


% colorbar

if colBar
    cb = colorbar;
    cb.Ticks = [0 1];
    cb.TickLabels = [ cLim ];
end


% anatomical outlines 

if drawLines
    % circumference
    hold on
    theta = linspace(0, 2*pi, 200);           % 200 points around the circle
    r     = coord_2d(nChan+1,2)/2;                   % radius based on nasion's coordinate
    cx    = 0.5; cy = 0.5;                       % center coordinates
    x     = cx + r*cos(theta);
    y     = cy + r*sin(theta);
    circumf = plot(x, y, 'LineWidth', lineWidth, 'Color',lineCol);
    uistack(circumf,'bottom')

    axis equal

    % add left and right ears
    leln = line([-0.03  -0.025],[0.38 0.56],'LineWidth', lineWidth, 'Color',lineCol);
    reln = line([ 1.03   1.025],[0.38 0.56],'LineWidth', lineWidth, 'Color',lineCol);
    % add nose
    lnln = line([0.4 0.5],[0.95 1]+0.05,'LineWidth', lineWidth, 'Color',lineCol);
    rnln = line([0.6 0.5],[0.95 1]+0.05,'LineWidth', lineWidth, 'Color',lineCol);

    uistack(leln,"bottom")
    uistack(reln,"bottom")
    uistack(lnln,"bottom")
    uistack(rnln,"bottom")

end




ax = gca;
ax.Color = [1 1 1 0];
ax.XAxis.Color = [1 1 1 0];
ax.YAxis.Color = [1 1 1 0];
ax.Colormap = colMap;
axis square; 
axis equal; 
%title(projectionType);

hold off