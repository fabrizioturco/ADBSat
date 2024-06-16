function h = plot_surfq_LVLH_withWind(fileIn, modIn, yaw_deg, pitch_deg, roll_deg, param, view_flag, view_az, view_el,yawIndex,pitchIndex,rollIndex)
% Plots the surface mesh with color proportional to the chosen parameter
%
% Inputs:
%       file_name   : Name of the file containing the results (fiName_eqmodel)
%       modIn       : Path to modesl
%       yaw_deg     : Yaw angle [deg]
%       pitch_deg   : Pitch angle [deg]
%       roll_deg    : Roll angle [deg]
%       param       : Chosen visualization parameter
%       view_flag   : Flag for View
%       view_az     : Flag for azimuth
%       view_el     : Flag for elevation
%       yawIndex    : Index for local yaw angle in assembled matrix
%       pitchIndex  : Index for local pitch angle in assembled matrix
%       rollIndex   : Index for local roll angle in assembled matrix
%
% Outputs:
%       h           : A patch object, containing the shape colour-coded by the chosen parameter
%
% Author: David Mostaza-Prieto
% The University of Manchester
% December 2012
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru and Sabrina Livadiotti
%
% This file is part of the ADBSat toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------

% Load model mesh
[~,modName,~] = fileparts(modIn);
load(modIn);
x_g = meshdata.XData;
y_g = meshdata.YData;
z_g = meshdata.ZData;

% Load results for indicated yaw, pitch, and roll
s = load(fileIn);

% Convert to Radians
yaw     = deg2rad(yaw_deg);
pitch   = deg2rad(pitch_deg);
roll    = deg2rad(roll_deg);

% Coordinate transformation matrices
L_gb = [1 0 0; 0 -1 0; 0 0 -1]; % Body to Geometric
L_fb = [-1 0 0; 0 1 0; 0 0 -1]; % Body to Flight
L_LVLHw = eye(3);
L_bLVLH = angle2dcm(yaw, pitch, roll); % LVLH to Body
L_gw    = L_gb * L_bLVLH * L_LVLHw;    % Wind to Geometric
L_gLVLH = L_gb * L_bLVLH;

axlength = 1*s.aedb.LenRef;

hFig = figure;
hold on
axis equal
grid on

color_param = s.aedb.(param);
coloring = zeros(1,size(color_param,2));    
for i=1:size(color_param,2) % iterates over panels
    coloring(1,i) = color_param{1,i}(pitchIndex,yawIndex,rollIndex);
end

h = patch(x_g, y_g, z_g, coloring);
colorbar
string1 = strcat(param,' Surface Distribution');

string2 = strcat('Yaw: ',mat2str(yaw_deg),' deg, Pitch: ', mat2str(pitch_deg), ' deg, Roll: ', mat2str(roll_deg));
xlabel('x_{LVLH}'); ylabel('y_{LVLH}'); zlabel('z_{LVLH}')
title(char(string1,string2))


% Geometric to body system
rotate(h, [1 0 0], 180)
set ( gca, 'ydir', 'reverse' )
set ( gca, 'zdir', 'reverse' )

% Body to LVLH system
rotate(h, [1 0 0], roll*180/pi);
rotate(h, [0 1 0], pitch*180/pi);
rotate(h, [0 0 1], yaw*180/pi);

if view_flag==1
    view(view_az,view_el)
end

%v_rel
vdirX = s.aedb.vdir_X(pitchIndex,yawIndex,rollIndex);
vdirY = s.aedb.vdir_Y(pitchIndex,yawIndex,rollIndex);
vdirZ = s.aedb.vdir_Z(pitchIndex,yawIndex,rollIndex);
vdir_LVLH = inv(L_gLVLH) * [vdirX;vdirY;vdirZ];
v_rel   = quiver3(axlength,0,0,0.5*axlength*vdir_LVLH(1),0.5*axlength*vdir_LVLH(2),0.5*axlength*vdir_LVLH(3),color='r',LineWidth=1,MaxHeadSize=1);
l       = legend(v_rel,'$\mathbf{v}_{rel}$','Location','NorthWest', 'interpreter', 'latex');
set(l,'FontSize',13);

axis equal
xlim([-axlength,axlength])
ylim([-axlength,axlength])
zlim([-axlength,axlength])

% Set custom update function
dcm = datacursormode(hFig);
set(dcm,'UpdateFcn',{@myupdatefcn,coloring,param});
end

function txt = myupdatefcn(~,evt,data,name)
pos = get(evt,'Position');
ind = ceil(get(evt, 'DataIndex')/3);
txt = { sprintf('(x,y,z): (%g, %g, %g)', pos(1:3)),...
    sprintf('index: %g', ind),...
    sprintf('%s value: %g', name, data(ind))
    };
end

%------------- END OF CODE --------------
