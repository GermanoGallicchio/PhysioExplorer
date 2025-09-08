function pe_y1x2View(pe_cfg, results, L, R, viewParams)
% ... Description ...
%
% INPUT:
%
%
% OUTPUT:
%
% Author:

%% shortcuts

ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;
ny1x2 = [ny1 nx2];

chanlocs = pe_cfg.dimensions.z3_chanLocs;

%% default parameters (if not entered)

fieldNames = fieldnames(viewParams);

% keep track of what is implicitly used by default and print it is pe_cfg.verbose = true
enteredByDefault = repmat("",0,0); % initialize

% colorMap
if ~any(strcmp(fieldNames,'colorMap'))
    enteredByDefault = [enteredByDefault "colorMap"];
    viewParams.colorMap = turbo;
end

% color limits
if ~any(strcmp(fieldNames,'cLim'))
    enteredByDefault = [enteredByDefault "cLim"];
    viewParams.cLim = prctile(results.statVal_obs(:),[1 99]);
end

% color bar visibility
if ~any(strcmp(fieldNames,'colorBarVisibility'))
    enteredByDefault = [enteredByDefault "colorBarVisibility"];
    viewParams.colorBarVisibility = true;
end

% color bar ticks
if viewParams.colorBarVisibility
    cTicks = unique([0 viewParams.cLim],'sorted');
    viewParams.cTicks = cTicks;
end




% contour mask
if ~any(strcmp(fieldNames,'contourMask'))
    enteredByDefault = [enteredByDefault "contourMask"];
    viewParams.contourMask = false;
end

% contour parameters
if viewParams.contourMask
    
    % contour color
    if ~any(strcmp(fieldNames,'contourColor'))
        enteredByDefault = [enteredByDefault "contourColor"];
        %viewParams.contourColor = [0 0 0];
        viewParams.contourColor = viewParams.colorMap(round(prctile(1:size(viewParams.colorMap,1),[75 25])),:);
    end

    % contour line width
    if ~any(strcmp(fieldNames,'contourLineWidth'))
        enteredByDefault = [enteredByDefault "contourLineWidth"];
        viewParams.contourLineWidth = 1.5;
    end
end




% x axis

if ~any(strcmp(fieldNames,'hLim'))
    enteredByDefault = [enteredByDefault "hLim"];
    viewParams.hLim = prctile(pe_cfg.dimensions.x2_vec,[0 100]);
end

if ~any(strcmp(fieldNames,'xAxisVisibility'))
    enteredByDefault = [enteredByDefault "xAxisVisibility"];
    viewParams.xAxisVisibility = 'true';
end

if ~any(strcmp(fieldNames,'xLabel'))
    enteredByDefault = [enteredByDefault "xLabel"];
    viewParams.xLabel = [pe_cfg.dimensions.x2_lbl ' ' '[' pe_cfg.dimensions.x2_units ']'];
end






% y axis

if ~any(strcmp(fieldNames,'vLim'))
    enteredByDefault = [enteredByDefault "vLim"];
    viewParams.vLim = prctile(pe_cfg.dimensions.y1_vec,[0 100]);
end

if ~any(strcmp(fieldNames,'yScale'))
    enteredByDefault = [enteredByDefault "yScale"];
    viewParams.yScale = 'linear';
end

if ~any(strcmp(fieldNames,'yAxisVisibility'))
    enteredByDefault = [enteredByDefault "yAxisVisibility"];
    viewParams.yAxisVisibility = 'true';
end

if ~any(strcmp(fieldNames,'yLabel'))
    enteredByDefault = [enteredByDefault "yLabel"];
    viewParams.yLabel = [pe_cfg.dimensions.y1_lbl ' ' '[' pe_cfg.dimensions.y1_units ']'];
end





if pe_cfg.verbose  &&  ~isempty(enteredByDefault)
    disp('Note: the following fields of viewParams were added by default')
    disp(enteredByDefault)
    disp('Check below the "viewParams" used by the code to see what fields can be changed')
    disp(viewParams)
end

%% highlight area mask
if viewParams.contourMask
    switch [pe_cfg.objective ' & ' pe_cfg.analysis]
        case 'permutationH0testing & empiricalL1_FDR'
            mask(:,:,:,1) = reshape(results.resampling.pVal_emp_FDR<pe_cfg.p_crit & results.statVal_obs>0, ny1, nx2, nz3); % positive
            mask(:,:,:,2) = reshape(results.resampling.pVal_emp_FDR<pe_cfg.p_crit & results.statVal_obs<0, ny1, nx2, nz3); % negative

        case 'permutationH0testing & theoreticalL1_clusterMaxT'
            mask(:,:,:,1) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs>0, ny1, nx2, nz3); % positive
            mask(:,:,:,2) = reshape(results.clusters.clusterMembership_mass_obs_Corrected~=0 & results.statVal_obs<0, ny1, nx2, nz3); % negative


            % --new section-- experimental does not work
            %         statVal_obsMasked = reshape(results.statVal_obs, ny1, nx2, nz3);
            %         masked(:,:,:,1) = statVal_obsMasked .* ~mask(:,:,:,1);
            %         masked(:,:,:,2) = statVal_obsMasked .* ~mask(:,:,:,2);
            %         mask = masked;
            % ---

        case 'bootstrapStability & empiricalL1_FDR'
            mask(:,:,:,1) = reshape(results.resampling.inference.BR_rob>2, ny1, nx2, nz3); % positive
            mask(:,:,:,2) = reshape(results.resampling.inference.BR_rob<-2, ny1, nx2, nz3); % negative

        case 'bootstrapStability & theoreticalL1_clusterMaxT'
            error('not yet coded')

        otherwise
            error('not yet coded')
    end
end

%% decide values to plot
% depending on results content 

xVals = pe_cfg.dimensions.x2_vec;
yVals = pe_cfg.dimensions.y1_vec;

% TO DO: split into two options: results based or L based
values2plot = reshape(results.statVal_obs, ny1, nx2, nz3);

%% figure

% draw
im = imagesc(xVals,yVals,values2plot);

% colormap
im.Parent.Colormap = viewParams.colorMap;

% color limits
im.Parent.CLim = viewParams.cLim;


% x-axis
im.Parent.XLim = viewParams.hLim;
if any(strcmp(fieldNames,'xTick'))
    xTick = viewParams.xTick;
    im.Parent.XTick = xTick;
end
if viewParams.xAxisVisibility
    im.Parent.XAxis.Label.String = viewParams.xLabel;
else
    im.Parent.XAxis.TickLabels = [];
end
im.Parent.XMinorTick = "off";

% y-axis
im.Parent.YLim = viewParams.vLim;
if any(strcmp(fieldNames,'yTick'))
    yTick = viewParams.yTick;
    im.Parent.YTick = yTick;
end
im.Parent.YDir = 'normal';
im.Parent.YScale = viewParams.yScale;
if viewParams.yAxisVisibility
    im.Parent.YAxis.Label.String = viewParams.yLabel;
else
    im.Parent.YAxis.TickLabels = [];
end
im.Parent.YMinorTick = "off";

% title
if any(strcmp(fieldNames,'title'))
    titleTxt = viewParams.title;
    title(titleTxt);
end



% contour
if viewParams.contourMask
    hold on
    for mIdx = 1:size(mask,4)
        if nnz(mask(:,:,:,mIdx)) > 0
            contour2plot = sum(mask(:,:,:,mIdx),4);
            col = viewParams.contourColor(mIdx,:);
            col = col.^(1/4);
            [~,ct] = contour(xVals,yVals,contour2plot,1, ...
                'LineWidth',viewParams.contourLineWidth);
            ct.EdgeColor = col;
            ct.ShowText = 'off';
            
        end
    end
end

% colorbar
if viewParams.colorBarVisibility 
    cb = colorbar;
    cb.Limits = viewParams.cLim;
    cb.Ticks = viewParams.cTicks;
    
    if any(strcmp(fieldNames,'colorBarLabel'))
        cbTxt = viewParams.colorBarLabel;
        cb.Label.String = cbTxt;
    end
    
end





