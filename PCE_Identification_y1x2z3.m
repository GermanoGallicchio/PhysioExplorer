function [clusterMatrix, clusterSize, clusterMass, clusterRobMass] = ...
    ClusterAnalysis_Identification_y1x2z3(statMatrix, pvalMatrix, DimStruct, p_crit, pixelSign, figFlag)

% This function identifies clusters in three dimensions: y1, x2, z3.
% The nature of the dimensions is described in the DimStruct structure.
% Two examples are below.
%   example #1: y1=Frequency, x2=Time, z3=Channel
%   example #2: y1=Frequency4Amplitude, x2=Frequency4Phase, z3=Channel
% 
%
% INPUT: 
%
% statMatrix:       matrix of statistical values (e.g., rho, tval) 
%                   dimensions: 1=chan, 2=freq, 3=time
% 
% pvalMatrix        matrix of pvalues associated with statMatrix 
%                   dimensions: 1=chan, 2=freq, 3=time
%
% neighborMatrix    matrix (chan x chan format)
%                   1=neighbors, 0=no
%
% DimStruct         structure defining the dimensions. Following example #1 above
%                     DimStruct.y1_lbl      = 'Freq';  % label for dimension y1
%                     DimStruct.y1_contFlag = 1;       % the variable is on one continuous scale 1=yes, 0=no
%                     DimStruct.y1_vec      = freqVec; % vector of values. needed only if contFlat==1
%                     DimStruct.y1_units    = 'Hz';    % label for the units of dimension y1. needed only if contFlat==1
%                     DimStruct.x2_lbl      = 'Time';   % label for dimension x2
%                     DimStruct.x2_contFlag = 1;        % the variable is on one continuous scale 1=yes, 0=no
%                     DimStruct.x2_vec      = timeVec;  % vector of values. needed only if contFlat==1
%                     DimStruct.x2_units    = 's';      % label for the units of dimension x2. needed only if contFlat==1
%                     DimStruct.z3_lbl      = 'Channel';
%                     DimStruct.z3_contFlag = 0;
%                     DimStruct.z3_chanlocs = EEG.chanlocs; % channel locations structure in the eeglab format
%                     DimStruct.z3_neighborMatrix = neighborMatrix; created by ClusterAnalysis_ChannelNeighborhood.m
%
% p_crit            1 digit represeting critical p value (e.g., 0.05) for clustering purposes only
%
% pixelSign         1 or -1 to search, respectively, for clusters with positive or negative statitical values
%
% figFlag           1=plot figure, 0=don't
%
% OUTPUT:
%
% clusterMatrix 
%
% clusterSize 
% 
% clusterMass 
% 
% clusterRobMass
%
% 
% Lay-term exaplanation: 
% Conceptually, imagine a 3d space made of "cubes" or "voxels". Two
% significant voxels belong in the same cluster if they share a whole
% "face". In other words, if they have the same coordinates for two
% dimensions and the other dimension is contiguous or neighboring.
%
% Acknowledgments: 
% inspired by Groppe's find_clusters.m and its subfunction follow_clust.m
% https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox 
% 
% written by Germano Gallicchio 
% germano.gallicchio@gmail.com

%% determine the nature of the dimensions

% determine the analysis design (types of variable) 
dimensions      = ["y1" "x2" "z3"];
presentInDesign = zeros(size(dimensions));
for dimIdx = 1:length(dimensions)
    presentInDesign(dimIdx) = logical(isfield(DimStruct, dimensions(dimIdx)+"_lbl"));
end
Label = repmat("",size(dimensions));
isContinuous = zeros(size(dimensions));
isChannel = zeros(size(dimensions));
for dimIdx = 1:length(dimensions)
    if presentInDesign(dimIdx)
        Label(dimIdx) = DimStruct.(dimensions(dimIdx)+"_lbl");
        isContinuous(dimIdx) = DimStruct.(dimensions(dimIdx)+"_contFlag");
        if isContinuous(dimIdx)==0
            isChannel(dimIdx) = contains(DimStruct.(dimensions(dimIdx)+"_lbl"),'channel','IgnoreCase',true);
        end
    end
end


% create dimension table
dimensions = dimensions(:);
presentInDesign = presentInDesign(:);
Label = Label(:);
isContinuous = isContinuous(:);
isChannel = isChannel(:);
DimTable = [table(dimensions)  table(Label)  table(presentInDesign) table(isContinuous) table(isChannel)];
% give feedback to user
% disp(DimTable(DimTable.presentInDesign==1,:))
% fprintf("dimensions to use: " +join(DimTable.dimensions(DimTable.presentInDesign==1)) + "\n")

designLabel = join(DimTable.Label(DimTable.presentInDesign==1),'_');

%% sanity checks

if DimStruct.z3_contFlag==0
    %disp('z3 is not continuous (e.g., EEG channels)')
    chanlocs = DimStruct.z3_chanlocs;
    neighborMatrix = DimStruct.z3_neighborMatrix;
else
    error('z3 is continuous (e.g., time, freq). great but not yet done coded')
end

% check that stat matrix has the expected dimensions
dimVec_obs = size(statMatrix);  % numerosity of each dimension (observed in statMatrix)
dimVec_exp = [length(DimStruct.y1_vec)  length(DimStruct.x2_vec)  length(DimStruct.z3_chanlocs)]; % numerosity of each dimension (expected from DimStruct)
if contains(designLabel,'channel','IgnoreCase',true)
    dimVec_exp = dimVec_exp([3 1 2]);    % move channel to first dimension
end
if ~isequal(dimVec_exp, dimVec_obs)
    error('dimensionality incongruence between statMatrix and what expected based on DimStruct')
else
    dimVec = dimVec_obs;
    clear dimVec_obs dimVec_exp
end
switch designLabel
    case {"Freq_Time_Channel" "FreqEEG_FreqHP_Channel"}
        nz3 = dimVec(1);
        ny1 = dimVec(2);
        nx2 = dimVec(3);
    case {"Freq_Time" "Freq_Freq" "Time_Time"}
        ny1 = dimVec(1);
        nx2 = dimVec(2);
end


% check that stat and pval matrices have the same dimensionality
if ~isequal(size(statMatrix), size(pvalMatrix))
    error('statMatrix and pvalMatrix must have the same size')
end

% check that neighborMatrix is square
if contains(designLabel,'channel','IgnoreCase',true)
    if ~isequal(size(neighborMatrix,1), size(neighborMatrix,2))
        error('neighborMatrix must be a square matrix')
    end
end

% check that statMatrix first dimension is channels
if contains(designLabel,'channel','IgnoreCase',true)
    if ~isequal(size(statMatrix,1), size(neighborMatrix,2))
        error('size(statMatrix,1) must be same as size(neighborMatrix,1)')
    end
end

if abs(pixelSign)~=1
    error('clusterSign must be 1 (positive) or -1 (negative)')
end

if ~(figFlag==0 | figFlag==1)
    error('figFlag must be 1=draw figure or 0=don''t')
end

%% get data

% initialize clusterMatrix
clusterMatrix = zeros(size(statMatrix));

%%

% apply threshold (based on p < p_crit) and separate pos from neg pixels
threshMatrix = pvalMatrix<p_crit;
threshMatrix = threshMatrix .* sign(statMatrix); % distinguish positive and negative clusters by their sign (+ or -)

% find idx of Above-Threshold Pixels (positive or negative, depending on pixelSign)
ATPixels = find(threshMatrix==pixelSign);

% get ND grids corresponding with the above-treshold pixels
switch designLabel
    case {"Freq_Time_Channel"  "FreqEEG_FreqHP_Channel"}
        [z3Matrix, y1Matrix, x2Matrix] = ndgrid(1:nz3,1:ny1,1:nx2);
    case {"Freq_Time" "Time_Time" "Freq_Freq"}
        [y1Matrix, x2Matrix] = ndgrid(1:ny1,1:nx2);
end

%% implementation

% from here below, "voxels" is used for "above-threshold voxel"
nClust = 0; % initialize / clusters numerosity
voxel_registrationStatus  = false(1,length(ATPixels)); % 0=not belonging to any cluster, 1=member of a cluster
voxel_registrationCluster = nan(1,length(ATPixels));   % number corresponding to the cluster they belong to

for pixel_idx = 1:length(ATPixels)    % loop through above-threshold pixels

    % if this cluster is not registered, register it in the next available cluster
    if voxel_registrationStatus(pixel_idx)==false
        nClust = nClust + 1;
        voxel_registrationStatus(pixel_idx)  = true;
        voxel_registrationCluster(pixel_idx) = nClust;
    end

    % search for neighbors
    toggleSwitch = true; % keep searching until this is true
    while toggleSwitch==true
        recruitCounter = 0;

        for recruiter_idx = find(voxel_registrationCluster==nClust)
            for recruitable_idx = find(voxel_registrationStatus==0)

                % y1 (eg, freq) distance between recruiter and recruitable
                dist_y1 = abs(y1Matrix(ATPixels(recruiter_idx)) - y1Matrix(ATPixels(recruitable_idx)));


                % x2 (eg, time) distance between recruiter and recruitable
                dist_x2 = abs(x2Matrix(ATPixels(recruiter_idx)) - x2Matrix(ATPixels(recruitable_idx)));
                
                % z3 (eg, channel) "distance" between recruiter and recruitable
                if z3Matrix(ATPixels(recruiter_idx))==z3Matrix(ATPixels(recruitable_idx)) % they are from the same channel
                    dist_z3 = 0;
                else
                    if neighborMatrix(z3Matrix(ATPixels(recruiter_idx)),z3Matrix(ATPixels(recruitable_idx)))==1  % they are from neighboring channels
                        dist_z3 = 1;
                    else % they are not from neighboring channels
                        dist_z3 = 2;
                    end
                end
                    
                % compound distance
                dist_compound = dist_y1 + dist_x2 + dist_z3;

                % decide on this voxel recruitment outcome
                % to join the cluster, at least two of three dist vars need to be 0 and the other 1
                if dist_compound<=1  % the voxel joins this cluster
                    voxel_registrationCluster(recruitable_idx) = nClust;
                    voxel_registrationStatus(recruitable_idx) = true;
                    recruitCounter = recruitCounter + 1;
                end
            end
        end
        
        % continue as long as there is a new recruit
        if recruitCounter==0
            toggleSwitch = false;
        end

    end

end

% update clusterMatrix
clusterMatrix(ATPixels) = voxel_registrationCluster;
clusterMatrix = clusterMatrix*pixelSign;  % if clusters are negative, make the registration number negative too

%% compute cluster statistics (e.g., size, mass, robust mass)

clusterSize    = nan(1,nClust); % initialize
clusterMass    = nan(1,nClust); % initialize
clusterRobMass = nan(1,nClust); % initialize
for clIdx = 1:nClust
    idx = clusterMatrix==clIdx*pixelSign; % idx of pixels belonging to this cluster

    clusterSize(clIdx)    = sum(idx(:));
    clusterMass(clIdx)    = mean(statMatrix(idx))   * clusterSize(clIdx); % can be positive or negative depending on content (ie, interest in pos or neg clusters)
    clusterRobMass(clIdx) = median(statMatrix(idx)) * clusterSize(clIdx); % can be positive or negative depending on content (ie, interest in pos or neg clusters)

    % sanity check: all stat values within this cluster have the same sign
    if pixelSign*sum(sign(statMatrix(idx)))~=sum(idx(:))
        error('this cluster contains a mix of positive and negative stats values. likely a bug in the code')
    end

end

%% post-processing sanity check
% plot

if figFlag==1
    
    % channel x freq x time plot (ie, many 2-d plots)
    figure(6); clf
    lyt = tiledlayout('flow');
    for chanIdx = 1:nz3
        nexttile(lyt)
        mat2plot = squeeze(clusterMatrix(chanIdx,:,:));
        imagesc(1:nx2,1:ny1,mat2plot)
        title(['chan ' num2str(chanIdx)])
    end
    sgtitle(['num of clusters: ' num2str(nClust)])


    % scalp maps
    figure(7); clf
    y1Idx = randi(ny1);  % pick a time point at random
    x2Idx = randi(nx2);  % pick a time point at random
    radius = 0.70; % for plotting reasons only (radius is not in the computations)
    nexttile()
    surrogateData = randn(1,nz3);
    % all chans
    topoplot(surrogateData, chanlocs, 'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0.5, 'emarker', {'o','k',3,1}); hold on;  % plot all channels
    clusters = nonzeros(unique(clusterMatrix(:,y1Idx,x2Idx)));
    tmpCol = lines(length(clusters));
    counter = 0;
    for clusterIdx = clusters'
        counter = counter +1;
        % cluster chans
        chan2plot = find(clusterMatrix(:,y1Idx,x2Idx)==clusterIdx);
        topoplot(surrogateData(chan2plot), chanlocs(chan2plot),  'electrodes', 'on', 'style', 'blank', 'plotrad',radius, 'headrad', 0, 'emarker', {'.',tmpCol(counter,:) ,15,5}); hold on; % plot neighbors of seed channels
    end
    title({['y1-dim point: ' num2str(y1Idx)]  ['x2-dim point: ' num2str(x2Idx)]})
    colormap(lines)
end
