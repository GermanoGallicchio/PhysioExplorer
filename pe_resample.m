function rowIdx = pe_resample(pe_cfg)

% Descr
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

nIterations = pe_cfg.nIterations;
nRow = height(pe_cfg.designTbl);

%% extract info 

% num of unique subjects
[unique_subjID,~,subjIdxPerRow] = unique(pe_cfg.designTbl.subjID,'stable');
nUnique_subjID = length(unique_subjID);

% num of subjID repeats
nRepeats_subjID = nRow / nUnique_subjID;

% unique groups
unique_groupID = string(unique(pe_cfg.designTbl.groupID,'stable'));
nUnique_groupID = length(unique_groupID);

% num of repeated-measure levels factors
rmF_col_idx = contains(pe_cfg.designTbl.Properties.VariableNames,'rmFactor');
nUnique_RMFactors = nnz(rmF_col_idx);
nUnique_RMLevelsPerFactors = nan(1,nUnique_RMFactors);
colCounter = 0;
for rmIdx = find(rmF_col_idx)
    colCounter = colCounter + 1;
    nUnique_RMLevelsPerFactors(1,colCounter) = length(unique(pe_cfg.designTbl{:,rmIdx}));
end
nUnique_RMLevels = prod(nUnique_RMLevelsPerFactors);

designTbl = pe_cfg.designTbl;

%% map subjects to groups and to rows

% group belonging for each unique subject
groupBelonging = strings(nUnique_subjID,1);
for subjIdx = 1:nUnique_subjID
    mask       = designTbl.subjID == unique_subjID(subjIdx);
    groupBelonging(subjIdx) = designTbl.groupID(find(mask,1));
end

% rows per each subject
rowIdxBySubj = accumarray(subjIdxPerRow, (1:nRow)', [], @(x){x});

% original group row-lists and counts
groupRowsByG  = cell(nUnique_groupID,1);
origGroupRowCnt = zeros(nUnique_groupID,1);
for gIdx = 1:nUnique_groupID
    subs = find(groupBelonging==unique_groupID(gIdx));
    groupRowsByG{gIdx}    = vertcat(rowIdxBySubj{subs});
    origGroupRowCnt(gIdx) = numel(groupRowsByG{gIdx});
end


%% build indices for resampling rows of L and/or R
% note: the last iteration is always the observed order

rng(pe_cfg.randomSeed); % random seed generator

switch [pe_cfg.objective ' & ' num2str(pe_cfg.designCode)]

    case {'permutationH0testing & 1  0  0'  
            'permutationH0testing & 0  1  0'
            'permutationH0testing & 1  1  0'}
        % correlation, group comparison, group comparison with covariate

        % H0 testing: build row permutation pattern for hypothesis(es), depending on the analysis design
        
        rowIdx = zeros(nRow,nIterations);
        for itIdx = 1:nIterations
            if itIdx==1
                rowIdx(:,itIdx) = 1:nRow; % same order for first iteration
            else
                rowIdx(:,itIdx) = randperm(nRow)'; % permuted order
            end
        end

    case {'permutationH0testing & 0  0  1'}
        % condition comparison
        
        condperm = @(x) (x(randperm(nRepeats_subjID))); % initialize condition swapping function

        rowIdx = zeros(nRow,nIterations);  % initialize
        for itIdx = 1:nIterations
            if itIdx==1
                rowIdxBySubj_perm = rowIdxBySubj; % same order for the first iteration
            else
                rowIdxBySubj_perm = cellfun(condperm,rowIdxBySubj,'UniformOutput',false); % permute order within each subject
            end
            for s = 1:nUnique_subjID
                rows = rowIdxBySubj{s}; % rows associated with this subject
                rowIdx(rows,itIdx) = rowIdxBySubj_perm{s}; % now associates permuted rows
            end
        end
    
        % sanity check: each subject occupies the same rows s/he was occupiying initially
        for itIdx = 1:nIterations
            for s = 1:nUnique_subjID
                rows_orig = rowIdxBySubj{s};
                rows_perm = rowIdx(rows_orig,itIdx);
                if ~length(intersect(rows_orig,rows_perm))==nUnique_RMLevels
                    error('likely coding bug: permuted sample has different subject to row mapping than original')
                end
            end
        end
        
    
    case {  'bootstrapStability & 0  0  1'
            'bootstrapStability & 0  1  0'
            'bootstrapStability & 0  1  1'
            'bootstrapStability & 1  0  0'
            'bootstrapStability & 1  0  1'
            'bootstrapStability & 1  1  0'
            'bootstrapStability & 1  1  1'
            }
        % ..description..
        rowIdx = zeros(nRow,nIterations);  % initialize

        for itIdx = 1:nIterations

            if itIdx==1
                % original order
                rowIdx(:, itIdx) = 1:size(rowIdx,1);
            else
                rowCounter = 1;  % row counter
                for gIdx = 1:nUnique_groupID
                    subsInG = find(groupBelonging==unique_groupID(gIdx));
                    nSubs   = length(subsInG);

                    % draw subjects with replacement (bootstrap)
                    draws = subsInG(randi(nSubs, [nSubs,1]));

                    for d = 1:nSubs
                        rows = rowIdxBySubj{draws(d)}; % rows associated with this subjects
                        sz   = numel(rows); % num of rows for this subject

                        rowIdx(rowCounter:(rowCounter-1+sz), itIdx) = rows;

                        rowCounter  = rowCounter + sz;
                    end
                end
            end
        end
        
        % sanity check 1: group counts preserved in bootstrap samples
        for itIdx = 1:nIterations
            idx = rowIdx(:,itIdx);
            for gIdx = 1:nUnique_groupID
                sampCnt = sum(ismember(idx, groupRowsByG{gIdx}));
                if ~sampCnt==origGroupRowCnt(gIdx)
                    error('likely coding bug: bootstrap sample has different group numerosity(ies) than original')
                end
            end
        end

        % sanity check 2: no subject carries an unexpected number of rows
        for itIdx = 1:nIterations
            idx = rowIdx(:,itIdx);
            blockSizes   = cellfun(@numel, rowIdxBySubj); % how many rows each subject would bring if they were sampled
            subjectCounts = cellfun(@(rows) sum(ismember(idx, rows)), rowIdxBySubj); % how many rows of each subject actually appear in the bootstrap sample
            partial = rem(subjectCounts, blockSizes) ~= 0; % flag any subject whose count is not an integer multiple of its blockSize
            if any(partial)
                error('likely coding bug: at least one subject did not carry all their conditions in the bootstrap sample')
            end
        end

    otherwise
        error(['not yet coded: ' pe_cfg.objective ' ' num2str(pe_cfg.designCode)])

end


          