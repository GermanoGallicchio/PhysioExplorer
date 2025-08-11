function [U,V] = pe_signConvention(Lz,Rz,U,V)

% INPUT:
%   Lz           - 
%
%   Rz           - 
%
%   U
%   V
%
%
% OUTPUT:
%
%
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)


% correlation sign between latent variables and dominant variable per mode (mIdx)
for mIdx = 1:size(V,2)

    % L-block latent variable for mode mIdx
    tL = Lz * U(:,mIdx);        % n×1

    % dominant variable for mode mIdx
    [~, domIdx] = max(abs(U(:,mIdx))); % find which original L‐column dominated U(:,mIdx)
    L_dom = Lz(:, domIdx);   % n×1  % that original variable time‐series

    % correlate the two vectors
    % if negative, flip both U and V
    if corr(tL, L_dom) < 0
        U(:,mIdx) = -U(:,mIdx);
        V(:,mIdx) = -V(:,mIdx);

        %fprintf(' just out of curiosity: sign change ')
    end
end