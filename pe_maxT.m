function results = pe_maxT(pe_cfg,results)

% Descr
%
% INPUT:
%
%   pe_cfg
%
%
% OUTPUT:
%
% results
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

nIterations = pe_cfg.nIterations;
metrics = results.resampling.metrics;
p_crit = pe_cfg.p_crit;

% num of clusters or modes
metrics_lbl = fieldnames(results.resampling.metrics);
metrics_lbl = metrics_lbl(~strcmp(metrics_lbl,'id')); % remove the 'id' column if present (it's present for cluster analysis)
metrics_num = length(metrics_lbl);
switch pe_cfg.analysis
    case 'theoreticalL1_clusterMaxT'
        nClust = length(metrics(end).(metrics_lbl{1}));
        r = nClust;
    case 'PLS_SVD'
        nModes = results.lvl2.nModes;
        r = nModes;
end

%% maxT correction

% most extreme value within each iteration (maxT)
metrics_maxT = metrics;  % initialize
for itIdx = 1:nIterations

    for metIdx = 1:metrics_num
        % all values within this iteration
        metricVal = [metrics(1,itIdx).(metrics_lbl{metIdx})]; 

        % find the largest
        [maxVal, maxIdx] = max(abs(metricVal));
        if maxVal==0; keyboard; end

        if ~isempty(maxVal)
            % metrics: keep only the largest 
            metrics_maxT(1,itIdx).(metrics_lbl{metIdx}) = metricVal(maxIdx);

            % (if present) 'id': keep only the largest 
            if any(strcmp(fieldnames(metrics_maxT),'id'))  &&  metIdx==1
                metrics_maxT(1,itIdx).id = metrics_maxT(1,itIdx).id(maxIdx);
            end
        else
            
            % metrics: set it zero
            metrics_maxT(1,itIdx).(metrics_lbl{metIdx}) = 0;

            % (if present and not done already) 'id': set it to none
            if any(strcmp(fieldnames(metrics_maxT),'id'))  &&  metIdx==1
                metrics_maxT(1,itIdx).id = NaN;
            end
        end
        
    end
end

% compute threshold corresponding with one tail p_crit
for metIdx = 1:metrics_num
    varLbl = metrics_lbl{metIdx};
    inference_maxT.thresholds.(varLbl) = prctile(abs([metrics_maxT.(metrics_lbl{metIdx})]),100-p_crit *100);
end


% compute p values
for msIdx = 1:metrics_num
    H0distribution = [metrics_maxT.(metrics_lbl{msIdx})];
    obsVal = metrics(1).(metrics_lbl{msIdx}); % observed measure (e.g., cluster mass, singular value)
    pval_maxT = nan(1,length(obsVal));
    for nIdx = 1:length(obsVal)
        pval_maxT(1,nIdx) = sum(abs(H0distribution)>=abs(obsVal(nIdx))) / length(H0distribution);
    end
    inference_maxT.pval.(metrics_lbl{msIdx}) = pval_maxT;
end

% view the distribution
if pe_cfg.figFlag  &&  nIterations>2
    figure(); clf;
    fig = gcf; fig.Units = "normalized"; fig.Position = [0.05 0.2 0.9 0.6];
    tld = tiledlayout('flow');
    yLbl = ['num of simulated samples (out of ' num2str(nIterations) ' iterations)'];
    tld.YLabel.String = yLbl;
    tld.Title.String = 'maxT null distribution';
    for metIdx = 1:metrics_num
        nexttile(tld);

        % H0 distribution histogram
        histogram(abs([metrics_maxT.(metrics_lbl{metIdx})]), 'FaceColor', [1 0 0],'FaceAlpha',0.5);
        xlabel(metrics_lbl{metIdx})

        % H0 distribution significance threshold
        hold on
        xline(inference_maxT.thresholds.([metrics_lbl{metIdx}]), 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 2);

        % observed data
        xline([metrics(1).(metrics_lbl{metIdx})],'LineStyle','-.','Color',[0.5 0.5 0]);

        % H0 distribution (all points)
        yyaxis right
        nVals = length([metrics_maxT.(metrics_lbl{metIdx})]);
        plot([metrics_maxT.(metrics_lbl{metIdx})],randn(1,nVals)+ones(1,nVals),'x','Color',[0 0 1]);
        set(gca,'YTick',[])
        legend({metrics_lbl{metIdx} [metrics_lbl{metIdx} '_{threshold}'] [metrics_lbl{metIdx} '_{observed}']},'Location','Best')

    end
end

% add to results
switch pe_cfg.analysis 
    case 'theoreticalL1_clusterMaxT'
        results.clusters.inference_maxT = inference_maxT;
        results.clusters.metrics_maxT = metrics_maxT;
    case 'PLS_SVD'
        results.pls_svd.inference_maxT = inference_maxT;
        results.pls_svd.metrics_maxT = metrics_maxT;
end