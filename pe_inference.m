function results = pe_inference(pe_cfg,results)

% Descr
% permutation (hypothesis-testing inference) => pvalues (1 tail)
% bootstrap (estimation-based inference) => BR, 2-tail 95%CI (2 tails)
%
% empiricalL1_FDR
%   permutation - feature-wise inference // multiple-comparison corrected: FDR over selected dimensions
%   bootstrap   - feature-wise inference // multiple-comparison uncorrected 
%
% cluster
%   permutation - cluster-level inference // multiple-comparison corrected: maxT
%   bootstrap   - cluster-level inference (mass or size) // multiple-comparison uncorrected (but fewer comparisons than feature-wise) 
%
% pls_svd
%   permutation - mode-level inference // multiple-comparison corrected: maxT and perMode
%   bootstrap   - loading-level inference // multiple comparison uncorrected
%
% INPUT:
%
%   pe_cfg
%
%
% OUTPUT:
%
% rowIdx
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% shortcuts

p_crit = pe_cfg.p_crit;
nIterations = pe_cfg.nIterations;

ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;

R_ignore = pe_cfg.R_ignore;

%% implementation

switch [pe_cfg.objective ' & ' pe_cfg.analysis]
    case 'permutationH0testing & empiricalL1_FDR'

        % find how "extreme" each statVal is compared with its distribution of statVals under H0
        pval_emp = nan(1,size(results.statVal_obs,2)); % initialize
        for colIdx = 1:size(results.statVal_obs,2)
            if pe_cfg.R_ignore(colIdx)
                continue
            end
            vals_perm = results.resampling.statVal_resamp(:,colIdx);
            val_obs = results.statVal_obs(colIdx);
            pval_emp(1,colIdx) = sum(abs(vals_perm)>=abs(val_obs)) / nIterations;
        end

        % --- TO DO: put in its own function ---
        % FDR correction
        pval_emp_FDR = nan(1,size(results.statVal_obs,2)); % initialize
        FDRdim_idx = find(pe_cfg.FDR_params.dimensions); % dimensions to pool for the correction
        otherdim_idx = setdiff(1:3, FDRdim_idx, 'stable'); % the other dimensions

        % permute to bring upfront the dimensions over which pvalues will be pooled
        perm  = [FDRdim_idx, otherdim_idx];
        pval_emp_permuted = permute(reshape(pval_emp,[ny1 nx2 nz3]),perm);
        R_ignore_permuted = permute(reshape(R_ignore,[ny1 nx2 nz3]),perm);

        % reshape to matrix to have columns of p values upon which the correction is applied (one iteration per column)
        sz = size(pval_emp_permuted);
        L  = prod(sz(1:numel(FDRdim_idx)));
        M  = prod(sz(numel(FDRdim_idx)+1:end));
        pval_emp_matrix = reshape(pval_emp_permuted, L, M);
        R_ignore_matrix = reshape(R_ignore_permuted, L, M);

        % apply FDR-BH along each column
        pval_emp_FDR_matrix = nan(size(pval_emp_matrix));
        for colIdx = 1:M
            pvalVec    = pval_emp_matrix(:,colIdx);
            RignoreVec = R_ignore_matrix(:,colIdx);
            pvalVec2use = pvalVec(~RignoreVec); % remove the p values corresponding with features that are not of interest in this analysis
            pvalVec2use_FDR = mafdr(pvalVec2use , 'BHFDR', true);  % FDR BH correction
            pval_emp_FDR_matrix(~RignoreVec,colIdx) = pvalVec2use_FDR; % collocate pvalues in the longer matrix where they belong
        end

        % reshape back
        pval_emp_FDR_permuted = reshape(pval_emp_FDR_matrix, sz);
        pval_emp_FDR = reshape(ipermute(pval_emp_FDR_permuted, perm),[1 ny1*nx2*nz3]);

        % figure
        if pe_cfg.figFlag
            figure(); clf
            f = gcf; f.Units = 'normalized'; f.Position = [0.2    0    0.4    0.9];
            vertJitterVec = randn(1,length(pval_emp(~pe_cfg.R_ignore)));
            lp1 = semilogx(pval_emp(~pe_cfg.R_ignore),vertJitterVec.*ones(1,length(pval_emp(~pe_cfg.R_ignore))),'o');
            lp1.Parent.XLim = [0 1];
            lp1.Parent.XTick = [0 0.01 0.05 0.1 0.2 0.5 1];
            lp1.Parent.XAxis.Label.String = 'p value';
            lp1.Parent.XMinorTick = 'off';
            lp1.Parent.YAxis.Visible = 'off';
            hold on
            lp2 = semilogx(pval_emp_FDR(~pe_cfg.R_ignore),vertJitterVec.*ones(1,length(pval_emp(~pe_cfg.R_ignore))),'x');
            ln = xline(pe_cfg.p_crit);
            legend([lp1 lp2 ln],["uncorrected" "FDR corrected" "p_{crit}"],'Location','southoutside')
            title('pvalues before vs after FDR correction')
        end

        % add to results
        results.resampling.pVal_emp     = pval_emp;
        results.resampling.pVal_emp_FDR = pval_emp_FDR;
        % --- ---

    case 'permutationH0testing & theoreticalL1_clusterMaxT'
        results = pe_maxT(pe_cfg,results); % maxT

    case 'permutationH0testing & PLS_SVD'
        results = pe_maxT(pe_cfg,results); % maxT

        % --- mode-specific then FDR (only for pls_svd)
        % TO DO: make its own function like pe_maxT
        % compute threshold corresponding with one tail p_crit

        % initialize
        inference_perMode = struct(); % probably can be deleted once I created the new function for pe_??? p vale per each mode

        for metIdx = 1:metrics_num

            if size([metrics(1).(metrics_lbl{metIdx})],2)~=r
                continue % this metric is not recorded once per mode
            end
            keyboard
            varLbl = metrics_lbl{metIdx};
            allVals = [metrics.(metrics_lbl{metIdx})];

            for rIdx = 1:r
                inference_perMode.thresholds.([varLbl '_mode' num2str(rIdx)]) = prctile(abs(allVals(rIdx,:)),100-p_crit *100);
            end
        end

        % compute p values
        for msIdx = 1:metrics_num
            H0distribution = [metrics.(metrics_lbl{msIdx})];
            if size(H0distribution,1)~=r
                continue % this metric is not recorded once per mode
            end
            pval_perMode = nan(1,r);
            for rIdx = 1:r
                obsVal = metrics(1).(metrics_lbl{msIdx})(rIdx); % observed measure (e.g., singular value)
                pval_perMode(1,rIdx) = sum(abs(H0distribution(rIdx,:))>=abs(obsVal)) / size(H0distribution,2);
            end
            inference_perMode.pval.(metrics_lbl{msIdx}) = pval_perMode;
            % FDR correction (Benjamini–Hochberg algorithm)
            pval_perMode_FDR = mafdr(pval_perMode, 'BHFDR', true);
            inference_perMode.pval_FDR.(metrics_lbl{msIdx}) = pval_perMode_FDR;
        end
        % add to results
        results.lvl2.inference_perMode = inference_perMode;
        % --- ---


        
    case 'bootstrapStability & empiricalL1_FDR'
        % compute bootstrap inference
        BR     = mean(results.statVal_obs,1) ./ std(results.resampling.statVal_resamp,0,1);
        BR_rob = mean(results.statVal_obs,1) ./ mad(results.resampling.statVal_resamp,1,1)*1.4826;
        CIlo = prctile(results.resampling.statVal_resamp(:,:),(pe_cfg.p_crit*100)/2,1);      % lower bound
        CIup = prctile(results.resampling.statVal_resamp(:,:),100-(pe_cfg.p_crit*100)/2,1);  % upper bound
        % apply R_ignore
        BR(logical(pe_cfg.R_ignore))     = NaN;
        BR_rob(logical(pe_cfg.R_ignore)) = NaN;
        CIlo(logical(pe_cfg.R_ignore))   = NaN;
        CIup(logical(pe_cfg.R_ignore))   = NaN;
        % add to results
        results.resampling.inference.BR     = BR;
        results.resampling.inference.BR_rob = BR_rob;
        results.resampling.inference.CIlo   = CIlo;
        results.resampling.inference.CIup   = CIup;
    case 'bootstrapStability & theoreticalL1_clusterMaxT'
        
        keyboard % TO DO: use observed value as numerator for BR and BR_rob // probably it will be by the time you are writing this: results.lvl1.statVal_obs
        BR     = mean(results.statVal_obs,1) ./ std(results.resampling.statVal,0,1);
        BR_rob = mean(results.statVal_obs,1) ./ mad(results.resampling.statVal,1,1)*1.4826;
        CIlo = prctile(results.lvl1.statVal(:,:),(pe_cfg.p_crit*100)/2,1);      % lower bound
        CIup = prctile(results.lvl1.statVal(:,:),100-(pe_cfg.p_crit*100)/2,1);  % upper bound
        % add to results
        results.lvl1.inference.BR     = BR;
        results.lvl1.inference.BR_rob = BR_rob;
        results.lvl1.inference.CIlo   = CIlo;
        results.lvl1.inference.CIup   = CIup;

    case 'bootstrapStability & pls_svd'

        % level 1 (loading per mode wise)
        % R-block scores
        BR = results.lvl1.V_obs(:,:) ./ std(results.lvl1.V_boot(:,:,:),0,3);
        BR(logical(pe_cfg.R_ignore),:) = 0; % apply R_ignore
        BR_rob = results.lvl1.V_obs(:,:) ./ mad(results.lvl1.V_boot(:,:,:),1,3) * 1.4826;
        BR_rob(logical(pe_cfg.R_ignore),:) = 0; % apply R_ignore
        CIlo = prctile(results.lvl1.V_boot(:,:,:),(pe_cfg.p_crit*100)/2,3);      % lower bound
        CIup = prctile(results.lvl1.V_boot(:,:,:),100-(pe_cfg.p_crit*100)/2,3);  % upper bound
        % add to results
        results.lvl1.inference.BR        = BR;
        results.lvl1.inference.BR_robust = BR_rob;
        results.lvl1.inference.CIlo      = CIlo;
        results.lvl1.inference.CIup      = CIup;
        % does it make sense to do the same for the L-block scores?
end

%% tidy up
% after resampling inference is made, we don't need to keep the resampling
% iterations in memory

switch [pe_cfg.objective ' & ' pe_cfg.analysis]
    case 'permutationH0testing & empiricalL1_FDR'
        results.resampling.statVal_resamp = 'deleted2saveMemory';

    case 'permutationH0testing & theoreticalL1_clusterMaxT'
        results = pe_clusterPruning(pe_cfg,results);
        results.resampling.metrics = 'deleted2saveMemory';

    case 'permutationH0testing & PLS_SVD'
        results.lvl2.metrics = results.lvl2.metrics(1,end);

    case 'bootstrapStability & empiricalL1_FDR'
        results.resampling.statVal_resamp = 'deleted2saveMemory';

    case 'bootstrapStability & theoreticalL1_clusterMaxT'
        results.lvl1.statVal = results.lvl1.statVal(end,:);
        results.lvl1.pVal    = results.lvl1.pVal(end,:);

    case 'bootstrapStability & PLS_SVD'            
                
        results.lvl1.U_boot = 'deleted2saveMemory';
        results.lvl1.V_boot = 'deleted2saveMemory';
        results.lvl1.LU_boot = 'deleted2saveMemory';
        results.lvl1.RV_boot = 'deleted2saveMemory';
end




