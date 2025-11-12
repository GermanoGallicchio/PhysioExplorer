function rowIdx = pe_resample(pe_cfg)
% SYNTAX:
%      rowIdx = pe_resample(pe_cfg)
%
% DE
% Build resampling/permutation indices for rows in designTbl based on the
% objective and design code. Supports permutation H0 testing and bootstrap
% stability while preserving grouping and within-subject structure.
%
% INPUT (fields of pe_cfg):
%   - designTbl: table with variables:
%       subjID  : subject identifier (scalar per row)
%       groupID : group identifier (scalar per row)
%       rmFactor* (optional): repeated-measure factors (one column per factor)
%   - nIterations: positive integer, number of resamples; column 1 is original order
%   - objective : 'permutationH0testing' | 'bootstrapStability'
%   - designCode: 1x3 numeric vector (e.g., [1 0 0], [0 1 0], [0 0 1])
%   - randomSeed (optional): integer seed for reproducibility
%
% OUTPUT:
%   - rowIdx: nRow x nIterations matrix of row indices (double)
%
% Notes:
%   - For 'permutationH0testing & 0  0  1' (condition comparison), rows are
%     permuted within subjects to preserve the within-subject block.
%   - For bootstrap, subjects are sampled with replacement within each group
%     and all of their rows (conditions) are carried over together.
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% basic validation and shortcuts
if ~isfield(pe_cfg, 'designTbl') || ~istable(pe_cfg.designTbl)
    error('pe_cfg.designTbl must be a table');
end
designTbl = pe_cfg.designTbl;

if ~all(ismember({'subjID','groupID'}, designTbl.Properties.VariableNames))
    error('designTbl must contain subjID and groupID columns');
end

if ~isfield(pe_cfg,'nIterations') || ~isnumeric(pe_cfg.nIterations) || ~isscalar(pe_cfg.nIterations) || pe_cfg.nIterations<1
    error('pe_cfg.nIterations must be a positive scalar');
end

if ~isfield(pe_cfg,'objective') || (~ischar(pe_cfg.objective) && ~isstring(pe_cfg.objective))
    error('pe_cfg.objective must be a string/char');
end

if ~isfield(pe_cfg,'designCode') || ~isnumeric(pe_cfg.designCode) || numel(pe_cfg.designCode)~=3
    error('pe_cfg.designCode must be a 1x3 numeric vector');
end

nIterations = double(pe_cfg.nIterations);
nRow = height(designTbl);

% set RNG if provided
if isfield(pe_cfg,'randomSeed') && ~isempty(pe_cfg.randomSeed)
    rng(pe_cfg.randomSeed, 'twister');
end

%% extract info
% subjects
[unique_subjID,~,subjIdxPerRow] = unique(designTbl.subjID,'stable');
nUnique_subjID = numel(unique_subjID);

% groups
unique_groupID = string(unique(designTbl.groupID,'stable'));
nUnique_groupID = numel(unique_groupID);

% repeated measures
rmF_col_idx = contains(designTbl.Properties.VariableNames,'rmFactor');
nUnique_RMFactors = nnz(rmF_col_idx);
nUnique_RMLevelsPerFactors = nan(1,nUnique_RMFactors);
colCounter = 0;
for rmIdx = find(rmF_col_idx)
    colCounter = colCounter + 1;
    nUnique_RMLevelsPerFactors(1,colCounter) = numel(unique(designTbl{:,rmIdx}));
end
if isempty(nUnique_RMLevelsPerFactors)
    nUnique_RMLevels = 1; % no RM factors
else
    nUnique_RMLevels = prod(nUnique_RMLevelsPerFactors);
end

%% map subjects to groups and to rows
% group belonging for each unique subject
groupBelonging = strings(nUnique_subjID,1);
for subjIdx = 1:nUnique_subjID
    mask       = designTbl.subjID == unique_subjID(subjIdx);
    groupBelonging(subjIdx) = string(designTbl.groupID(find(mask,1)));
end

% rows per each subject (cell array)
rowIdxBySubj = accumarray(subjIdxPerRow, (1:nRow)', [], @(x){x});

% original group row-lists and counts (for bootstrap checks)
groupRowsByG    = cell(nUnique_groupID,1);
origGroupRowCnt = zeros(nUnique_groupID,1);
for gIdx = 1:nUnique_groupID
    subs = find(groupBelonging==unique_groupID(gIdx));
    groupRowsByG{gIdx}    = vertcat(rowIdxBySubj{subs});
    origGroupRowCnt(gIdx) = numel(groupRowsByG{gIdx});
end

%% build indices for resampling rows of L and/or R
% column 1 is always the observed (original) order
key = [char(pe_cfg.objective) ' & ' num2str(pe_cfg.designCode)];

switch key

    case {'permutationH0testing & 1  0  0', ...   % correlation
          'permutationH0testing & 0  1  0', ...   % group comparison
          'permutationH0testing & 1  1  0'}       % group comparison with covariate
        rowIdx = zeros(nRow, nIterations);
        rowIdx(:,1) = (1:nRow)';
        for itIdx = 2:nIterations
            rowIdx(:,itIdx) = randperm(nRow)';
        end

    case {'permutationH0testing & 0  0  1'}       % condition comparison (within subject)
        % permute within each subject block to preserve within-subject structure
        rowIdx = zeros(nRow, nIterations);
        rowIdx(:,1) = (1:nRow)';
        for itIdx = 2:nIterations
            rowIdxBySubj_perm = cellfun(@(x) x(randperm(numel(x))), rowIdxBySubj, 'UniformOutput', false);
            for s = 1:nUnique_subjID
                rows = rowIdxBySubj{s};
                rowIdx(rows,itIdx) = rowIdxBySubj_perm{s};
            end
        end

        % sanity check: each subject keeps the same set of rows
        for itIdx = 1:nIterations
            for s = 1:nUnique_subjID
                rows_orig = rowIdxBySubj{s};
                rows_perm = rowIdx(rows_orig,itIdx);
                if numel(intersect(rows_orig, rows_perm)) ~= nUnique_RMLevels
                    error('likely coding bug: permuted sample changed subject-to-row mapping')
                end
            end
        end

    case { 'bootstrapStability & 0  0  1', ...
           'bootstrapStability & 0  1  0', ...
           'bootstrapStability & 0  1  1', ...
           'bootstrapStability & 1  0  0', ...
           'bootstrapStability & 1  0  1', ...
           'bootstrapStability & 1  1  0', ...
           'bootstrapStability & 1  1  1' }

        rowIdx = zeros(nRow, nIterations);
        rowIdx(:,1) = (1:nRow)'; % original order

        for itIdx = 2:nIterations
            rowCounter = 1;
            for gIdx = 1:nUnique_groupID
                subsInG = find(groupBelonging==unique_groupID(gIdx));
                nSubs   = numel(subsInG);

                % sample subjects with replacement within group
                draws = subsInG(randi(nSubs, [nSubs,1]));

                for d = 1:nSubs
                    rows = rowIdxBySubj{draws(d)}; % all rows for this subject
                    sz   = numel(rows);
                    rowIdx(rowCounter:(rowCounter-1+sz), itIdx) = rows;
                    rowCounter  = rowCounter + sz;
                end
            end
        end

        % sanity check 1: group counts preserved in bootstrap samples
        for itIdx = 1:nIterations
            idx = rowIdx(:,itIdx);
            for gIdx = 1:nUnique_groupID
                sampCnt = sum(ismember(idx, groupRowsByG{gIdx}));
                if sampCnt ~= origGroupRowCnt(gIdx)
                    error('likely coding bug: bootstrap sample has different group numerosity than original')
                end
            end
        end

        % sanity check 2: subjects contribute whole blocks only
        blockSizes     = cellfun(@numel, rowIdxBySubj);
        for itIdx = 1:nIterations
            idx = rowIdx(:,itIdx);
            subjectCounts = cellfun(@(rows) sum(ismember(idx, rows)), rowIdxBySubj);
            partial = rem(subjectCounts, blockSizes) ~= 0; % non-multiples indicate partial blocks
            if any(partial)
                error('likely coding bug: at least one subject did not carry all their conditions in the bootstrap sample')
            end
        end

    otherwise
        error('not yet coded: %s %s', char(pe_cfg.objective), num2str(pe_cfg.designCode));
end


          