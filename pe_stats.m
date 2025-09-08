function results = pe_stats(pe_cfg, L, R)

% This is the front-end wrapper function that reads the input, calls the
% "back-end" functions to do the jobs (e.g., understand what design the
% user wants, do sanity checks, do the analysis) and presents the output to
% the user. 
%
% [results, clusterMetrics, clustThreshold] = PE_Stats_y1x2z3(dataArray, PE_parameters, DimStruct)
%
% INPUT:
%   L           - 
%
%   R           - 
%
%   pe_cfg
%
%
% OUTPUT:
%
%
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% sanity checks

% L and R are matrices
if ~length(size(L))==2
    error('L must be a matrix')
end
if ~length(size(R))==2
    error('R must be a matrix')
end

% L and R have the same num of rows
if ~size(L,1)==size(R,1)
    error('L and R must have the same num of rows')
end

% designTbl must exist
fieldNeeded = 'designTbl';
if ~any(contains(fieldnames(pe_cfg),fieldNeeded))
    error(['pe_cfg needs field: ' fieldNeeded])
end

% L has the same num of rows as in pe_cfg.designTbl
if ~size(L,1)==size(pe_cfg.designTbl,1)
    error('L and designTbl must have the same num of rows')
end

% univariate analysis only does one comparison at the time
if ismember(pe_cfg.analysis, ["empiricalL1_FDR" "theoreticalL1_clusterMaxT"])
    if size(L,2)~=1
        error(['more than one matrix L columns (i.e., more than one comparison) not supported for ' num2str(pe_cfg.designCode) ' ' pe_cfg.analysis ])
    end
end

% check objective field matches one of the options
objectiveList = ["permutationH0testing" "bootstrapStability"]';
if ~ismember(pe_cfg.objective,objectiveList)
    disp(objectiveList)
    error('pe_cfg.objective must be one one of the above')
end

% check analysis field matches one of the options
analysisList = ["empiricalL1_FDR" "theoreticalL1_clusterMaxT" "pls_svd"]';
if ~ismember(pe_cfg.analysis,analysisList)
    disp(analysisList)
    error('pe_cfg.analysis must be one one of the above')
end

%% apply row ignore

% matrices L and R
L = L(~pe_cfg.row_ignore,:);
R = R(~pe_cfg.row_ignore,:);

% designTbl
pe_cfg.designTbl = pe_cfg.designTbl(~pe_cfg.row_ignore,:);

%% keep a copy of the original matrices

R_orig   = R;
L_orig   = L;

%% shortcuts

[m, pL] = size(L);
[~, pR] = size(R);

nIterations = pe_cfg.nIterations;
ny1 = pe_cfg.dimensions.y1_num;
nx2 = pe_cfg.dimensions.x2_num;
nz3 = pe_cfg.dimensions.z3_num;

%% defaults value

% if group not entered in designTable: put 1s all over
varLbl = pe_cfg.designTbl.Properties.VariableNames;
if ~any(contains(varLbl,'groupID'))
    warning('"groupID" column in pe_cfg.designTable not entered. I am assuming you have one group.')
    disp('note: to achieve the same and not see the warning above create a group variable and use all 1s)')
    pe_cfg.designTbl.groupID = ones(size(pe_cfg.designTbl,1),1);
end

% if no rmFactor entered: put 1s all over
varLbl = pe_cfg.designTbl.Properties.VariableNames;
if ~any(contains(varLbl,'rmFactor'))
    warning('"rmFactor1" column in pe_cfg.designTable not entered. I am assuming you have one rmFactor with one level. (note: use all 1s if there is only one rmLevel)')
    pe_cfg.designTbl.rmFactor1 = ones(size(pe_cfg.designTbl,1),1);
end

%% parse design
% understand what analysis the user wants to do


% designOptions = pe_designOptions;  % can delete this row
designCode = pe_parseDesign(pe_cfg,L);
pe_cfg.designCode = designCode;
results.designCode = designCode;

%% perform the analysis

% get resampling indices
rowIdx = pe_resample(pe_cfg);  

switch [pe_cfg.analysis ' & ' num2str(pe_cfg.designCode)]
    case 'empiricalL1_FDR & 1  0  0' % -- correlation --

        % choose analysis subtype
        list = ["Pearson" "Spearman" "Kendall" "cylindrical"];
        [idx,tf] = listdlg('ListString',list,'SelectionMode','single','ListSize',[160 100],'PromptString','choose correlation type');
        if tf==0; idx=1; warning('You did not choose correlation type. I chose for you: Pearson'); end
        corrType = list(idx);
        if strcmp(corrType,'cylindrical')
            error('not coded yet') % it will require its own function to keep things tidy
        end

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);

            % perform test
            [statVal, pVal] = corr(L, R, 'type', corrType); % test
  
            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
            end

            if itIdx==1
                statVal_resamp = zeros(nIterations,pR);
            end
            statVal_resamp(itIdx,:) = statVal;

            pe_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs        = statVal_obs;
        results.pVal_obs           = pVal_obs;
        results.resampling.statVal_resamp = statVal_resamp;

    case 'empiricalL1_FDR & 0  1  0' % -- independent sample t-test --
        for itIdx = 1:nIterations
            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);

            % perform test
            varType = 'unequal';  % equal | unequal (for info see doc ttest2)
            [~,p,~,stats] = ttest2(R(L==max(L),:),R(L==min(L),:),'Vartype',varType);
            statVal = stats.tstat;
            pVal = p;

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
            end

            if itIdx==1
                statVal_resamp = zeros(nIterations,pR);
            end
            statVal_resamp(itIdx,:) = statVal;

            pe_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs        = statVal_obs;
        results.pVal_obs           = pVal_obs;
        results.resampling.statVal_resamp = statVal_resamp;

    case 'empiricalL1_FDR & 0  0  1' % -- paired sample t-test --
        for itIdx = 1:nIterations

            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);
            
           

            % rows of conditions belonging to first and second halves
            cond_firstHalf = 1:length(L)/2;
            cond_secondHalf = 1+length(L)/2:length(L);
            % rows of conditions with the largest and lowest constrast code
            % ie, find which hald corresponds with the largest constrast code
            contCodes = [unique(L(cond_firstHalf)) unique(L(cond_secondHalf))];
            [~, sortIdx] = sort(contCodes,'descend');
            if diff(sortIdx)<0
                cond_maxHalf = cond_secondHalf;
                cond_minHalf = cond_firstHalf;
            else
                cond_maxHalf = cond_firstHalf;
                cond_minHalf = cond_secondHalf;
            end

            % perform test
            [~,p,~,stats] = ttest(R(cond_maxHalf,:),R(cond_minHalf,:));
            statVal = stats.tstat;
            pVal = p;

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
            end

            if itIdx==1
                statVal_resamp = zeros(nIterations,pR);
            end
            statVal_resamp(itIdx,:) = statVal;

            pe_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs        = statVal_obs;
        results.pVal_obs           = pVal_obs;
        results.resampling.statVal_resamp = statVal_resamp;

    case 'theoreticalL1_clusterMaxT & 1  0  0' % -- correlation --

        % initialize lvl2 (cluster) metrics
        resampling  = repmat(struct('id', [], 'size', [], 'mass', []), 1, nIterations);

        % choose analysis subtype
        list = ["Pearson" "Spearman" "Kendall" "cylindrical"];
        [idx,tf] = listdlg('ListString',list,'SelectionMode','single','ListSize',[160 100],'PromptString','choose correlation type');
        if tf==0; idx=1; warning('You did not choose correlation type. I chose for you: Pearson'); end
        corrType = list(idx);
        if strcmp(corrType,'cylindrical')
            error('not coded yet') % it will require its own function to keep things tidy
        end

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);

            % perform test
            [statVal, pVal] = corr(L, R, 'type', corrType); % test

            % form clusters
            [clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, statVal ,pVal);

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
                clusterMembership_obs = clusterMembership ;
                clustIDList_obs       = clustIDList ;
                metrics_obs           = metrics;
            end

            switch pe_cfg.objective
                case 'permutationH0testing'
                    resampling(1,itIdx) = metrics;

                case 'bootstrapStability'
                    if itIdx==1
                        statVal_boot = zeros(nIterations,ny1*nx2*nz3);
                    end
                    statVal_boot(itIdx,:) = statVal;
            end

            pe_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        results.statVal_obs                    = statVal_obs;
        results.pVal_obs                       = pVal_obs;
        results.clusters.clusterMembership_obs = clusterMembership_obs;
        results.clusters.clustIDList_obs       = clustIDList_obs;
        results.clusters.metrics_obs           = metrics_obs;
        
        switch pe_cfg.objective
            case 'permutationH0testing'
                results.resampling.metrics           = resampling;
            case 'bootstrapStability'
                results.resampling.statVal_boot      = statVal_boot(:,:);
        end

    case 'theoreticalL1_clusterMaxT & 0  1  0' % -- independent sample t-test --

        % initialize lvl2 (cluster) metrics
        resampling  = repmat(struct('id', [], 'size', [], 'mass', []), 1, nIterations);

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);

            % perform test
            varType = 'unequal';  % equal | unequal (for info see doc ttest2)
            [~,p,~,stats] = ttest2(R(L==max(L),:),R(L==min(L),:),'Vartype',varType);
            statVal = stats.tstat;
            pVal = p;

            % form clusters
            [clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, statVal ,pVal);

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
                clusterMembership_obs = clusterMembership ;
                clustIDList_obs       = clustIDList ;
                metrics_obs           = metrics;
            end

            switch pe_cfg.objective
                case 'permutationH0testing'
                    resampling(1,itIdx) = metrics;

                case 'bootstrapStability'
                    if itIdx==1
                        statVal_boot = zeros(nIterations,ny1*nx2*nz3);
                    end
                    statVal_boot(itIdx,:) = statVal;
            end

            pe_counter(itIdx,nIterations)  % iteration counter
        end
                
        % collate results
        results.statVal_obs                    = statVal_obs;
        results.pVal_obs                       = pVal_obs;
        results.clusters.clusterMembership_obs = clusterMembership_obs;
        results.clusters.clustIDList_obs       = clustIDList_obs;
        results.clusters.metrics_obs           = metrics_obs;

        switch pe_cfg.objective
            case 'permutationH0testing'
                results.resampling.metrics           = resampling;
            case 'bootstrapStability'
                results.resampling.statVal_boot      = statVal_boot(:,:);
        end

    case 'theoreticalL1_clusterMaxT & 0  0  1' % -- paired sample t-test --

        % initialize lvl2 (cluster) metrics
        resampling  = repmat(struct('id', [], 'size', [], 'mass', []), 1, nIterations);

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);

             % rows of conditions belonging to first and second halves
            cond_firstHalf = 1:length(L)/2;
            cond_secondHalf = 1+length(L)/2:length(L);
            % rows of conditions with the largest and lowest constrast code
            % ie, find which hald corresponds with the largest constrast code
            contCodes = [unique(L(cond_firstHalf)) unique(L(cond_secondHalf))];
            [~, sortIdx] = sort(contCodes,'descend');
            if diff(sortIdx)<0
                cond_maxHalf = cond_secondHalf;
                cond_minHalf = cond_firstHalf;
            else
                cond_maxHalf = cond_firstHalf;
                cond_minHalf = cond_secondHalf;
            end

            % perform test
            [~,p,~,stats] = ttest(R(cond_maxHalf,:),R(cond_minHalf,:));
            statVal = stats.tstat;
            pVal = p;

            % form clusters
            [clusterMembership, clustIDList, metrics] = pe_clusterForming(pe_cfg, statVal ,pVal);

            if itIdx==1
                statVal_obs           = statVal;
                pVal_obs              = pVal;
                clusterMembership_obs = clusterMembership ;
                clustIDList_obs       = clustIDList ;
                metrics_obs           = metrics;
            end

            switch pe_cfg.objective
                case 'permutationH0testing'
                    resampling(1,itIdx) = metrics;

                case 'bootstrapStability'
                    if itIdx==1
                        statVal_boot = zeros(nIterations,ny1*nx2*nz3);
                    end
                    statVal_boot(itIdx,:) = statVal;
            end

            pe_counter(itIdx,nIterations)  % iteration counter
        end
                
        % collate results
        results.statVal_obs                    = statVal_obs;
        results.pVal_obs                       = pVal_obs;
        results.clusters.clusterMembership_obs = clusterMembership_obs;
        results.clusters.clustIDList_obs       = clustIDList_obs;
        results.clusters.metrics_obs           = metrics_obs;

        switch pe_cfg.objective
            case 'permutationH0testing'
                results.resampling.metrics           = resampling;
            case 'bootstrapStability'
                results.resampling.statVal_boot      = statVal_boot(:,:);
        end

    case {'PLS_SVD & 1  0  0'   'PLS_SVD & 0  1  0'   'PLS_SVD & 1  1  0'}
        nModes = min(rank(L),rank(R));

        for itIdx = 1:nIterations

            % sort rows as appropriate
            [L,R] = pe_sortRows(pe_cfg,L_orig,R_orig,rowIdx,itIdx);

            % mean center columns of L and R
            Lz = normalize(L,'center');
            Rz = normalize(R,'center');

            % zscore columns of L (optional)
            if pe_cfg.pls_svdParams.zscoringVec(1)
                Lz = zscore(L);
            end
            % zscore columns of R (optional)
            if pe_cfg.pls_svdParams.zscoringVec(2)
                Rz = zscore(R);
            end

            % apply ignore mask
            Rz(:,logical(pe_cfg.R_ignore)) = 0;

            % scale by frobenius norm (optional)
            % DELETE THIS SECTION UNLESS JUSTIFIED
            %                     if logical(pe_cfg.pls_svdParams.froFlag)
            %                         Lz = Lz / norm(Lz,'fro');
            %                         Rz = Rz / norm(Rz,'fro');
            %                     end

            % -- temp code for PCA---
            %                     % PCA to get explained variance
            %                     [Rcoeff, score, ~, ~, explained] = pca(Rz);
            %                     % choose k to capture at least XX% variance
            %                     cumVar = cumsum(explained);
            %                     Rk = find(cumVar >= 90, 1, 'first');
            %                     disp(Rk)
            %                     Rz_pca = score(:,1:Rk);   % m×k
            % --- ---

            % cross-product
            %                     C = (Lz'*Rz_pca)/(m-1);  % for PCA
            C = (Lz'*Rz)/(m-1);

            % SVD
            [U,S,V] = svd(C,"econ");

            % --- temporary code for PCA
            %                     V = Rcoeff(:, 1:Rk)*V;
            % ---

            % sanity check: nModes as expected
            if nModes~=size(S,1)
                error('unexpected number of SVD modes')
            end

            % resolve sign uncertainty
            % match SVD outcome with observed data
            % only for last iteration
            if itIdx==1
                [U, V] = pe_signConvention(Lz,Rz,U,V);
            end

            % match bootstrap with original sample (the 1st one)
            % TO DO
            % 1. mode matching
            % 2. Procrustes rotation
            % 3. sign matching

            % metrics: level-1 (feature based)
            LU = Lz*U;
            RV = Rz*V;

            % metrics: level-2 (SV based)
            s = diag(S); % singular value of each mode
            p = s.^2 / sum(s.^2); % proportion of covariance explained by each mode
            r = diag(corr(LU, RV)); % Pearson correlation between left and right latent variables (i.e., scores = data * singular vectors) per mode
            % store the above as col vectors
            s = s(:);
            p = p(:);
            r = r(:);


            % store original-data values (observed)
            if itIdx==1
                U_obs = U;
                V_obs = V;
                LU_obs = LU;
                RV_obs = RV;
                s_obs = s';
                p_obs = p';
                r_obs = r';
            end

            % store output of each iteration as appropriate
            switch pe_cfg.objective
                case 'permutationH0testing'
                    % initialize lvl2 (mode wise) metrics
                    if itIdx==1
                        resampling  = repmat(struct('s', [], 'inertia', [], 'wilk', [], 'sequential', []), 1, nIterations);
                    end
                    % fill each field per iteration
                    resampling(:,itIdx).s = s;
                    resampling(:,itIdx).inertia = sum(s);
                    resampling(:,itIdx).wilk = flipud(cumsum(flipud(s(:).^2)));
                    resampling(:,itIdx).sequential = cumsum(s(:).^2);

                case 'bootstrapStability'
                    % initialize lvl1 (loading wise) metrics
                    % store singular and latent vectors (all iterations)
                    if itIdx==1
                        U_boot = zeros(size(L,2),nModes,nIterations);
                        V_boot = zeros(size(R,2),nModes,nIterations);
                        LU_boot = zeros(size(L,1),nModes,nIterations);
                        RV_boot = zeros(size(R,1),nModes,nIterations);
                    end
                    U_boot(:,:,itIdx) = U;
                    V_boot(:,:,itIdx) = V;
                    LU_boot(:,:,itIdx) = LU;
                    RV_boot(:,:,itIdx) = RV;
            end

            pe_counter(itIdx,nIterations)  % iteration counter
        end

        % collate results
        % observed data
        results.PLS_SVD.nModes = nModes;
        results.PLS_SVD.s_obs = s_obs;
        results.PLS_SVD.p_obs = p_obs;
        results.PLS_SVD.r_obs = r_obs;
        results.PLS_SVD.U_obs = U_obs;
        results.PLS_SVD.V_obs = V_obs;
        results.PLS_SVD.LU_obs = LU_obs;
        results.PLS_SVD.RV_obs = RV_obs;

        switch pe_cfg.objective
            case 'permutationH0testing'
                results.resampling     = resampling;
            case 'bootstrapStability'
                % all iterations at lvl1
                results.resampling.U_boot   = U_boot(:,:,:);
                results.resampling.V_boot   = V_boot(:,:,:);
                results.resampling.LU_boot  = LU_boot(:,:,:);
                results.resampling.RV_boot  = RV_boot(:,:,:);
        end

    otherwise % -- any other design --
        error(['not yet coded: ' pe_cfg.objective ' ' num2str(pe_cfg.designCode)])
end

%% compute inferential metrics
% compute p values for permutation testing
% compute BR and CI for bootstrap stability

results = pe_inference(pe_cfg,results);
       

%% cluster descriptive metrics

results = pe_describeClusters(pe_cfg,results);

