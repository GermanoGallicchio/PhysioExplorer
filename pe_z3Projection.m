function coord_2d = pe_z3Projection(coord_3d, projectionType)

% Description
% 
%   INPUT:
%
%       coord_3d
%
%       projectionType
%
%
%   OUTPUT:
%
%       coord_2d
%
%
%   Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

radius = 1; % hard coded

%% sanity checks

projectionOptions = ["orthographic" "azimuthalEquidistant" "azimuthalConformal"];
if ~ismember(projectionType,projectionOptions)
    disp(projectionOptions)
    error('projection option must be one of the above')
end
%% implementetation

az  = deg2rad([coord_3d.sph_theta]);  % azimuth
el  = deg2rad([coord_3d.sph_phi]);    % elevation

% colatitude: angular distance from the "north pole" (i.e., vertex, Cz)
colat = (pi/2) - el;

% compute coordinates
switch projectionType
    case 'orthographic'
        [x3, y3, z3] = sph2cart(az, el, radius);
        % simply ignore z
        xOrtho = x3;
        yOrtho = y3;
        xCoord = xOrtho;
        yCoord = yOrtho;

    case 'azimuthalEquidistant'
        % treat colat as the planar radius
        [xTopo, yTopo] = pol2cart(az, radius * colat);
        xCoord = xTopo;
        yCoord = yTopo;

    case 'azimuthalConformal'
        % use stereographic radius
        rho = tan(colat ./ 2);
        [xStereo, yStereo] = pol2cart(az, radius * rho);
        xCoord = xStereo;
        yCoord = yStereo;

end

coord_2d = [xCoord' yCoord'];




