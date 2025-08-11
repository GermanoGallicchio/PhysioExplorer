function designCode = pe_parseDesign(pe_cfg, L)
% deduces from the information what design the user wants to do 

% INPUT:
%   pe_cfg
%
%   L
%
%
% OUTPUT:
%
%

% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


%% sanity checks

designTbl_varLbl = pe_cfg.designTbl.Properties.VariableNames;

% design table must include one subjID column
if ~any(contains(designTbl_varLbl,'subjID'))
    error('design table must include one subjID column')
end

% design table must include one groupID column
if ~any(contains(designTbl_varLbl,'groupID'))
    error('design table must include one groupID column (note: use all 1s if there is only one group)')
end

% design table must include at least one rmFactor column
if ~any(contains(designTbl_varLbl,'rmFactor'))
    error('design table must include at least one "rmFactor" column (note: use all 1s if there is only one rmLevel)')
end

%% put all information on the table

nIterations = pe_cfg.nIterations;
nRow = height(pe_cfg.designTbl);

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


% --- DELETE THIS SECTION (UNNEEDED) --- 
% 
% % num of repeated-measure levels factors
% rmF_col_idx = contains(pe_cfg.designTbl.Properties.VariableNames,'rmFactor');
% nUnique_RMFactors = nnz(rmF_col_idx);
% nUnique_RMLevelsPerFactors = nan(1,nUnique_RMFactors);
% colCounter = 0;
% for rmIdx = find(rmF_col_idx)
%     colCounter = colCounter + 1;
%     nUnique_RMLevelsPerFactors(1,colCounter) = length(unique(pe_cfg.designTbl{:,rmIdx}));
% end
% nUnique_RMLevels = prod(nUnique_RMLevelsPerFactors);
% 
% % num of columns of L 
% nLCol = size(L,2);
% 
% nIterations = pe_cfg.nIterations;
% nRow = height(pe_cfg.designTbl);
% 
% % num of unique subjects
% [unique_subjID,~,subjIdxPerRow] = unique(pe_cfg.designTbl.subjID,'stable');
% nUnique_subjID = length(unique_subjID);
% 
% % num of subjID repeats
% nRepeats_subjID = nRow / nUnique_subjID;
% 
% % unique groups
% unique_groupID = string(unique(pe_cfg.designTbl.groupID,'stable'));
% nUnique_groupID = length(unique_groupID);
% 
% % num of repeated-measure levels factors
% rmF_col_idx = contains(pe_cfg.designTbl.Properties.VariableNames,'rmFactor');
% nUnique_RMFactors = nnz(rmF_col_idx);
% nUnique_RMLevelsPerFactors = nan(1,nUnique_RMFactors);
% colCounter = 0;
% for rmIdx = find(rmF_col_idx)
%     colCounter = colCounter + 1;
%     nUnique_RMLevelsPerFactors(1,colCounter) = length(unique(pe_cfg.designTbl{:,rmIdx}));
% end
% nUnique_RMLevels = prod(nUnique_RMLevelsPerFactors);
% 
% designTbl = pe_cfg.designTbl;
% --- 

%% sanity checks

% number of subjID repeats must match num of unique overall RM levels
if ~nRepeats_subjID==nUnique_RMLevels
    error('likely user design error: the number of subjID repeats must match the number of unique repeated-measure levels (overall across factors)')
end

%% deduct intended design

designOptions = pe_designOptions;


if nRepeats_subjID==1
    if nUnique_groupID<=1
        designCode = [1 0 0];
    else
        if size(L,2) < (nUnique_groupID-1)
            designCode = [0 1 0];
            warning('likely user design error: not enough columns in L')
        elseif size(L,2) == (nUnique_groupID-1)
            designCode = [0 1 0];
        else
            designCode = [1 1 0];
        end
    end
else
    if nUnique_groupID<=1
        if size(L,2) < nUnique_RMLevels-1
            designCode = [0 0 1];
            warning('likely user design error: not enough columns in L')
        elseif size(L,2) == nUnique_RMLevels-1
            designCode = [0 0 1];
        else
            designCode = [1 0 1];
        end
    else
        if size(L,2) < (nUnique_groupID-1)+(nUnique_RMLevels-1)+((nUnique_groupID-1)*(nUnique_RMLevels-1))
            designCode = [0 1 1];
            warning('likely user design error: not enough columns in L')
        elseif size(L,2) == (nUnique_groupID-1)+(nUnique_RMLevels-1)+((nUnique_groupID-1)*(nUnique_RMLevels-1))
            designCode = [0 1 1];
        else
            designCode = [1 1 1];
        end
    end
end

%% plot structure of L

if pe_cfg.figFlag
    figure()
    imagesc(normalize(L,'range',[-1 1]))
    title('structure of L')
    axis equal
    box off
    axis off
    colormap(parula)

end

%fprintf(['intended design: \n' design '\n'])

%% sanity checks on R
% TO DO ONE DAY but not essential. It's the user's responsibility to prepare the matrices (e.g., mean center appropriately, mainly for PLS-SVD)

%% 
if pe_cfg.verbose
    rowIdx = (designOptions{:,"codeDec"}==bin2dec(num2str(designCode)))';
    fprintf(['\n' 'Design num is: ' num2str(designOptions{rowIdx,"codeDec"}) ', ' designOptions{rowIdx,"lbl"}{1} '\n'])
end

