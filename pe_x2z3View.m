function pe_view(pe_cfg, results, L, R, viewParams)
% ... Description ...
%
% INPUT:
%
%
% OUTPUT:
%
% Author:

%% default parameters (if not entered)

fieldNames = fieldnames(viewParams);

if ~any(strcmp(fieldNames,'xLabel'))
    warning('viewParams.xLabel not entered, so I will guess')
    viewParams.xLabel = [pe_cfg.dimensions.x2_lbl ' ' '[' pe_cfg.dimensions.x2_units ']'];
end

if ~any(strcmp(fieldNames,'yLabel'))
    warning('viewParams.yLabel not entered, so I will guess')
    viewParams.yLabel = [pe_cfg.dimensions.y1_lbl ' ' '[' pe_cfg.dimensions.y1_units ']'];
end

if ~any(strcmp(fieldNames,'xAxis4z3Lbl'))
    warning('viewParams.xAxis4z3Lbl not entered, so I will guess')
    viewParams.xAxis4z3Lbl = ["P7" "P8"];
end

if ~any(strcmp(fieldNames,'yAxis4z3Lbl'))
    warning('viewParams.yAxis4z3Lbl not entered, so I will guess')
    viewParams.yAxis4z3Lbl = ["P7" "F7"];
end

if ~any(strcmp(fieldNames,'title4z3Lbl'))
    warning('viewParams.title4z3Lbl not entered, so I will guess')
    viewParams.title4z3Lbl = string({pe_cfg.dimensions.z3_chanLocs.labels});
end


if ~any(strcmp(fieldNames,'lineColorPalette'))
    warning('viewParams.lineColorPalette not entered, so I will guess')
    tmpLines = lines(5);
    viewParams.lineColorPalette = tmpLines([4 5],:); % lines/waveforms
end

if ~any(strcmp(fieldNames,'areaColorPalette'))
    warning('viewParams.areaColorPalette not entered, so I will guess')
    tmpLines = lines(5);
    viewParams.areaColorPalette = tmpLines([2 1],:); % highlight area
end


if ~any(strcmp(fieldNames,'xyLength'))
    warning('viewParams.xyLength not entered, so I will guess')
    viewParams.xyLength  = [0.10  0.10];
end

if ~any(strcmp(fieldNames,'xyPadding'))
    warning('viewParams.xyPadding not entered, so I will guess')
    viewParams.xyPadding = [0.00  0.00];
end

if ~any(strcmp(fieldNames,'projectionType'))
    warning('viewParams.projectionType not entered, so I will guess')
    viewParams.projectionType = 'azimuthalEquidistant';
end

if ~any(strcmp(fieldNames,'projectionWarping'))
    warning('viewParams.projectionWarping not entered, so I will guess')
    viewParams.projectionWarping = [1 1];
end

if ~any(strcmp(fieldNames,'hLim'))
    error('include viewParams.hLim')
end
if ~any(strcmp(fieldNames,'hStep'))
    error('include viewParams.hStep')
end
if ~any(strcmp(fieldNames,'vLim'))
    error('include viewParams.vLim')
end

%% shortcuts

ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;
ny1x2 = [ny1 nx2];



vLim = viewParams.vLim;
hLim = viewParams.hLim;
hStep = viewParams.hStep;


xLabel = viewParams.xLabel;
yLabel = viewParams.yLabel;

xAxis4z3Lbl = viewParams.xAxis4z3Lbl;
yAxis4z3Lbl = viewParams.yAxis4z3Lbl;
title4z3Lbl = viewParams.title4z3Lbl;

lineColorPalette = viewParams.lineColorPalette;
areaColorPalette = viewParams.areaColorPalette;


xyLength  = viewParams.xyLength;
xyPadding = viewParams.xyPadding;


projectionType = viewParams.projectionType; 
projectionWarping = viewParams.projectionWarping;
%% decide content of horizontal and vertical axes

switch num2str(ny1x2>1)
    case num2str([1   0])
        error('not coded yet / probably I will never enable this feature. simply swap dimensions y1 and x2 earlier and the job is done')

    case num2str([0  1])
        hAxisDim = 'x2';
        xVals = pe_cfg.dimensions.([hAxisDim '_vec']);

    case num2str([1 1])
        error('not coded yet')
        vAxisDim = 'y1';
        yVals = pe_cfg.dimensions.([vAxisDim '_vec']);
        hAxisDim = 'x2';
        xVals = pe_cfg.dimensions.([hAxisDim '_vec']);

    otherwise
        error('dimensions y1 and x2 are both null. perhaps you want to use z3Plot?')
end


%% more shortcuts





%% highlight area mask
switch [pe_cfg.objective ' & ' pe_cfg.analysis]
    case 'permutationH0testing & empiricalL1_FDR'
        mask(:,:,:,1) = reshape(results.resampling.pVal_emp_FDR<pe_cfg.p_crit & results.statVal_obs>0, ny1, nx2, nz3); % positive
        mask(:,:,:,2) = reshape(results.resampling.pVal_emp_FDR<pe_cfg.p_crit & results.statVal_obs<0, ny1, nx2, nz3); % negative

    case 'permutationH0testing & theoreticalL1_clusterMaxT'
        mask(:,:,:,1) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs>0, ny1, nx2, nz3); % positive
        mask(:,:,:,2) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs<0, ny1, nx2, nz3); % negative

    case 'bootstrapStability & empiricalL1_FDR'
        mask(:,:,:,1) = reshape(results.resampling.inference.BR_rob>2, ny1, nx2, nz3); % positive
        mask(:,:,:,2) = reshape(results.resampling.inference.BR_rob<-2, ny1, nx2, nz3); % negative

    case 'bootstrapStability & theoreticalL1_clusterMaxT'
        error('not yet coded')

    otherwise
        error('not yet coded')
end

%% axis (subplot) locations 

% shortcuts
chanlocs = pe_cfg.dimensions.z3_chanLocs;

% 3d --> 2d projections
coord_3d = chanlocs;
coord_2d = pe_z3Projection(coord_3d, projectionType);

% coordinate warping (nonlinear scaling) to fit the drawing area better
xSign = sign(coord_2d(:,1));
ySign = sign(coord_2d(:,2));
coord_2d(:,1) = xSign.*(abs(coord_2d(:,1)).^(projectionWarping(1))); % 0.85
coord_2d(:,2) = ySign.*(abs(coord_2d(:,2)).^(projectionWarping(2))); % 1.11

% create position vector
Position = coord_2d;
Position(:,1) = normalize(Position(:,1),'range',[0+xyPadding(1) 1-xyLength(1)-xyPadding(1)]);
Position(:,2) = normalize(Position(:,2),'range',[0+xyPadding(2) 1-xyLength(2)-xyPadding(2)]);
Position = [Position(:,:) repmat(xyLength(1),nz3,1) repmat(xyLength(2),nz3,1) ];

%% draw

for chanIdx = 1:length(chanlocs)

    % create axis set (subplot)
    ax(chanIdx) = axes();
    ax(chanIdx).Units = 'normalized';
    ax(chanIdx).Position = Position(chanIdx,:);
    
    % plot data series/waveform of data corresponding with largest L code
    rowIdx = L==max(L);
    val2plot_tmp = reshape(R(rowIdx,:),[sum(rowIdx) ny1 nx2 nz3]);
    val2plot = shiftdim(mean(val2plot_tmp(:,:,:,chanIdx),1),1);
    lp1 = plot(ax(chanIdx), xVals, val2plot); hold on
    lp1.LineWidth = 1;

    % plot data series/waveform of data corresponding with smallest L code
    rowIdx = L==min(L);
    val2plot_tmp = reshape(R(rowIdx,:),[sum(rowIdx) ny1 nx2 nz3]);
    val2plot = shiftdim(mean(val2plot_tmp(:,:,:,chanIdx),1),1);
    lp2 = plot(ax(chanIdx), xVals, val2plot);
    lp2.LineWidth = 0.8;

    % colors
    lp1.Color = viewParams.lineColorPalette(1,:);
    lp2.Color = viewParams.lineColorPalette(2,:);


    % x-axis
    ax(chanIdx).XLim = hLim;
    xTick = unique([fliplr(0:-hStep:min(xVals))  0:hStep:max(xVals)],'stable');
    ax(chanIdx).XTick = xTick;
    if ismember(chanlocs(chanIdx).labels,xAxis4z3Lbl)
        ax(chanIdx).XLabel.String = xLabel;
    else
        ax(chanIdx).XTickLabel = [];
    end

    % y-axis
    ax(chanIdx).YLim = vLim;
    if ismember(chanlocs(chanIdx).labels,yAxis4z3Lbl)
        ax(chanIdx).YLabel.String = yLabel;
    else
        ax(chanIdx).YTickLabel = [];
    end

    % title channel label
    if ismember(chanlocs(chanIdx).labels,title4z3Lbl)
        title(chanlocs(chanIdx).labels,'VerticalAlignment','bottom')
    end

    % significance area
    for signIdx = 1:size(mask,4)
        mask2plot = mask(1,:,chanIdx,signIdx);

        % semi-transparent areas
        d_sigIdx = diff([0 mask2plot 0]);  % highlight starting and ending points (i.e., moving from 0 to 1 gives "1" (starting), moving from 1 to 0 gives "-1" (ending), everything else is 0 and indicates no changes)
        area_startIdx = find(d_sigIdx==1);
        area_endIdx = find(d_sigIdx==-1)-1;
        step = mean(diff(xVals))/2;
        for islandIdx = 1:length(area_startIdx)
            a = xVals(area_startIdx(islandIdx));
            b = xVals(area_endIdx(islandIdx));
            % add step to allow an area instead of a line
            a = a - step;
            b = b + step;
            % draw
            ar(signIdx) = area(ax(chanIdx),[a b], range(vLim)*ones(1,2), vLim(1), ...
                'FaceColor' , areaColorPalette(signIdx,:), ...
                'FaceAlpha' , 0.25, ...        % 25% opacity
                'EdgeColor' , 'none');         % no edge line
            uistack(ar(signIdx),'bottom');
        end
    end

    % plot (statVal)
%     statValAxis = false;
%     if statValAxis 
%     yyaxis right
%     val2plot_tmp = reshape(results.statVal_obs,[ny1 nx2 nz3]);
% %     val2plot_tmp = reshape(results.resampling.inference.BR_rob,[ny1 nx2 nz3]);
%     val2plot = shiftdim(mean(val2plot_tmp(:,:,chanIdx),1),1);
%     lp3 = plot(ax(chanIdx), erptime, val2plot, 'LineWidth',0.8); hold on
%     % colors
%     %lp3.Color = lineColorPal(3,:);
%     lp3.Color = 0.5*[1 1 1];
%     lp3.Color(4) = 0.5;
%     lp3.LineWidth = 0.1;
%     lp3.LineStyle = '-';
%     ax(chanIdx).YAxis(2).Limits = statValLim;
%     ax(chanIdx).YAxis(2).Color = lp3.Color;
%     % y-axis (right)
%     ax(chanIdx).YAxis(2).TickValues = [statValLim(1) 0 statValLim(2)];
%     if ismember(chanlocs(chanIdx).labels,["P8" "F8"])
%         ax(chanIdx).YLabel.String = 'statistic score';
%     else
%         ax(chanIdx).YTickLabel = [];
%     end
%     end


    box off
    ax(chanIdx).Color = [0 0 0   0];
    ax(chanIdx).TickLength = [0.05 0.025];
    
end

% legend
ax2 = axes();
ax2.Position = [0 0 1 1];
ax2.Visible = 'off';

switch num2str(results.designCode)
    case num2str([1 0 0])
        warning('not coded yet')
        keyboard
    case num2str([0 1 0])
        lp1_lbl = unique(pe_cfg.designTbl.groupID(L==max(L)));
        lp2_lbl = unique(pe_cfg.designTbl.groupID(L==min(L)));
    case num2str([0 0 1])
        lp1_lbl = unique(pe_cfg.designTbl.rmFactor1(L==max(L)));
        lp2_lbl = unique(pe_cfg.designTbl.rmFactor1(L==min(L)));
    otherwise
        error('not coded yet')
end
legend(ax2,[lp1 lp2 ],[ lp1_lbl lp2_lbl ],'Position',[0.05 0.88 0.125 0.08],'Box','off')
