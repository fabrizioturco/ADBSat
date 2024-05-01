% Calculates local and global coefficients for the triangular mesh geometry.
% Produces one or several output .m files with the object characteristics.
%
% Inputs:
%       fiName   : Name of the .mat file containing the meshdata structure
%                   XData
%                   YData
%                   ZData
%                   MatID
%                   Areas
%                   SurfN
%                   BariC
%                   Lref
%       yawArray    : Yaw angle(s) [rad]
%       pitchArray  : Pitch angle(s) [rad]
%       rollArray   : Roll angle(s) [rad]
%       phi         : angle between flight direction and position vector [rad]
%       inc         : inclination of orbit [rad]
%       eqmodel     : String containing the name of the equation to be used to calculate coefficients
%       v_windmzu   : Wind velocity in Meridian-Zonal-Up reference frame [m/s]
%       magnv_sat   : Satellite velocity magnitude [m/s]
%       param_eq    : Parameters associated to "eqmodel"
%       flag_shad   : Flag to perform shadow analysis (1=perform, 0=no)
%       flag_sol    : Flag to calculate solar wind coefficients
%       pathin      : Path for loop runs (only for loop sims)
%
% Outputs:
%       fileOut  : Output file
%           yaw      : Yaw angle [rad]
%           pitch    : Pitch angle [rad]
%           roll     : Roll angle [rad]
%           tauDir   : Direction of the tangent vectors (3xN)
%           delta    : Angles between the flow and the faces [rad] (1xN)
%           cp       : Face pressure coefficient (1xN)
%           ctau     : Face shear coefficient (1xN)
%           cd       : Face drag coefficient (1xN)
%           cl       : Face lift coefficient (1xN)
%           Cf_w     : Body force coefficient in wind axes (3x1)
%           Cf_f     : Body force coefficient in flight axes (3x1)
%           Cf_b     : Body force coefficient in body axes (3x1)
%           Cf_LVLH  : Body force coefficient in LVLH axes (3x1)
%           C_D      : Body drag coefficient relative to AreaProj
%           Cm_B     : Body moment coefficient in body axes (3x1)
%           Aref     : Reference area used in calculations [m^2]
%           AreaProj : Projected area to the flow [m^2]
%           Lref     : Refrence length used for calcualtions [m]
%           param_eq : Structure containing input GSI parameters
%
% where N = number of faces.
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
% 
% CHANGELOG
% 13 Feb 2024: Now considering angles of yaw, pitch, and roll ((Euler angles) instead ofangle of attack and sideslip 
% (F. Turco)
%
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

function [fileOut] = calc_coeff_eulerAngles(fiName, respath, yawArray, pitchArray, rollArray, phi, inc, v_windmzu, magnv_sat, param_eq, flag_shad, flag_sol, del, verb)

[~,matName,~] = fileparts(fiName);

% Load mesh parameters
load(fiName,'meshdata');
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
areas = meshdata.Areas;
surfN = meshdata.SurfN;
barC = meshdata.BariC;
Lref = meshdata.Lref;
matID = meshdata.MatID;

indexYaw = length(yawArray);
indexPitch = length(pitchArray);
indexRoll = length(rollArray);


%Waitbar
if verb
    h = waitbar(0,'Please wait...','Name','Calculating coefficients...');
end

% Create output folder if required
if (indexYaw*indexPitch*indexRoll) > 1
    foldname = strcat(matName,'_',datestr(now,30),'_',num2str(randi(1000)));
    mkdir(fullfile(respath,filesep,foldname));
    pathsav = fullfile(respath,foldname);
    aedb = 1;
else
    pathsav = fullfile(respath);
    aedb = 0;
end

% Values to save in output
var_out = { 'yaw';'pitch';'roll';'tauDir';'delta';'cp';'ctau';'cd';'cl';...
            'Cf_w';'Cf_f';'Cf_b';'Cf_LVLH';'C_D';'Cm_B';...
            'Aref';'AreaProj';'Lref';'param_eq';'shadow'};
if flag_sol
    var_out = [var_out;{'Cf_s';'Cm_S'}];
end


for ii = 1:indexYaw
    yaw = yawArray(ii);

    for jj = 1:indexPitch
        pitch = pitchArray(jj);

        for kk = 1: indexRoll
            roll = rollArray(kk);
            
            %L_wb = dcmbody2wind(aoa, aos); % Body to Wind
            % L_wb = [cos(aos)*cos(aoa), sin(aos), sin(aoa)*cos(aos);...
            %     -sin(aos)*cos(aoa), cos(aos), -sin(aoa)*sin(aos);...
            %     -sin(aoa), 0, cos(aoa)]; % Body to Wind

            L_gb = [1 0 0; 0 -1 0; 0 0 -1]; % Body to Geometric
            L_fb = [-1 0 0; 0 1 0; 0 0 -1]; % Body to Flight
            L_bLVLH = angle2dcm(yaw, pitch, roll); % LVLH to Body
            L_gLVLH = L_gb * L_bLVLH;

            L_nedmzu = [-1 0 0; 0 1 0; 0 0 -1]; % MZU to NED
            if v_windmzu(1) == 0
                windAngle = deg2rad(90);
            else
                windAngle = atan(v_windmzu(2)/v_windmzu(1));
            end
            L_wmzu = angle2dcm(windAngle,0,0); % MZU to Wind (Rotation around z)
            L_LVLHned = angle2dcm(inc,0,pi/2-phi);  % NED to LVLH
            L_gned = L_gLVLH*L_LVLHned; % NED to Geometric
            L_gw = L_gLVLH * L_LVLHned * L_nedmzu * inv(L_wmzu); % Wind to Geometric

            % Debugging velocities
            v_windw = L_wmzu*v_windmzu; % Wind velocity in Wind (Alternative: v_windw = [-sqrt(dot(v_windmzu,v_windmzu)); 0; 0];)
            v_windned = L_nedmzu * v_windmzu;
            v_windLVLH = L_LVLHned * v_windned;

            % Transforming velocities
            v_windg = L_gLVLH * v_windLVLH; % Wind velocity in Geometric
            v_satLVLH = [-magnv_sat;0;0];   % Satellite velocity in LVLH
            v_satg = L_gLVLH * v_satLVLH;   % Satellite velocity in Geometric
    
            % Flow direction
            v_relg = v_satg + v_windg;
            vdir = v_relg/norm(v_relg);
            %v_rel = [-1;0;0]; % Relative velocity of atmosphere in LVLH system
            %vdir = L_gw * v_rel;
            %vdir = vdir/norm(vdir);
            
            % Surface normals and angles between faces and flow
            vMatrix = [vdir(1)*ones(1,length(surfN(1,:)));...
                vdir(2)*ones(1,length(surfN(1,:)));...
                vdir(3)*ones(1,length(surfN(1,:)))];
            
            % Angles between flow and normals
            delta = real(acos(dot(-vMatrix,surfN)));
            
            uD = vMatrix; % Unit drag vector
            uL = -cross(cross(uD,surfN),uD)./vecnorm(cross(cross(uD,surfN),uD)); % Unit lift vector for each panel
            % Undefined uL panels (plates normal to flow)
            col = find(all(isnan(uL),1));
            uL(:,col) = -surfN(:,col);
            % Negative dot product of unit drag and lift vectors with surface normal vector (March 2019)
            param_eq.gamma = dot(-uD,surfN);
            param_eq.ell = dot(-uL,surfN);
            
            % Local flat plate coefficients
            [cp, ctau, cd, cl] = mainCoeff(param_eq, delta, matID);
            
            if flag_sol
                [cn, cs] = coeff_solar(delta, param_eq);
            end
            
            % Backwards facing panels
            areaB = areas;
            areaB(delta*180/pi>90) = 0;
            
            shadow = zeros(size(areas));
            shadow(areaB == 0) = 1;
            
            % Shadow analysis
            if flag_shad
                [shadPan] = shadowAnaly(x, y, z, barC, delta, L_gw);
                
                cp(shadPan) = 0;
                ctau(shadPan) = 0;
                cd(shadPan) = 0;
                cl(shadPan) = 0;
                areaB(shadPan) = 0;
                
                cn(shadPan) = 0;
                cs(shadPan) = 0;            
                shadow(shadPan) = 0.5;
            end
      
            % Areas
            AreaProj = areaB*cos(delta)'; % Projected
            AreaT = sum(areas); % Total
            Aref = AreaT/2; % ADBSat Reference
            
            % Shear direction
            tauDir = cross(surfN,cross(vMatrix,surfN)); % direction of the shear coefficient
            ModT = sqrt(dot(tauDir,tauDir));
            tauDir = tauDir./[ModT;ModT;ModT]; % normalise
            
            TF = isnan(tauDir(1,:)); % Check any non defined tauDir (angles 90,180)
            tauDir(:,TF) = 0;
            
            % Force coefficents
            Cf_g    = 1/(Aref).*(tauDir*(ctau'.*areas') - surfN*(cp'.*areas'));
            Cf_b    = L_gb' * Cf_g;
            Cf_w    = L_gw' * Cf_g;
            Cf_f    = L_fb * L_gb' * Cf_g;
            Cf_LVLH = L_gLVLH' * Cf_g;

            % Drag coefficient
            C_D = -Cf_w(1)*Aref/AreaProj;
            
            % Moment coefficients
            Cmom_G = 1/(Aref*Lref).*(cross(barC,tauDir)*(ctau'.*areas') +...
                cross(barC,-surfN)*(cp'.*areas'));
            %Cmom_G_b = sum((1/(Aref*Lref).*(cross(barC, (tauDir.*(ctau.*areas) - surfN.*(cp.*areas))))),2);
            Cm_B = L_gb' * Cmom_G;
            
            % Solar coefficients
            if flag_sol
                cstau = cs.*sin(delta); % resolve incident to shear
                cs_n = cs.*cos(delta); % resolve incident to normal
                Cs_G = 1/(Aref).*(tauDir*(cstau'.*areas') - surfN*((cn+cs_n)'.*areas'));
                Cf_s = L_gw'*Cs_G;
                Cms_G = 1/(Aref*Lref).*(cross(barC,tauDir)*(cstau'.*areas') +...
                    cross(barC,-surfN)*((cn+cs_n)'.*areas'));
                Cm_S = L_gb'*Cms_G;
            end
            
            % Save file
            if aedb
                fileOut = fullfile(pathsav,[matName,'_yaw',mat2str(yaw*180/pi),...
                                                    '_pitch',mat2str(pitch*180/pi),...
                                                    '.roll',mat2str(roll*180/pi),...
                                                    '.mat']);
            else
                fileOut = [pathsav,'.mat'];
            end
            
            save(fileOut, var_out{:})
            
            if verb
                percent = (ii-1)*1/indexYaw + (jj-1)*1/(indexPitch*indexYaw) + kk*1/(indexPitch*indexYaw*indexRoll);
                waitbar(percent,h,sprintf('%6.4g',percent*100));
            end
        end   
    end
end

if verb
    delete(h); % clean-up waitbar
end

if (indexYaw*indexPitch*indexRoll) > 1
    % Creates a merged aerodynamic database from multiple .mat files
    fileOut  = mergeAEDB(pathsav, matName, del);
end


% Debugging reference frames
debugRF = 1;
if debugRF
    Fig = figure;
    hold on
    
    axlength = 2;

    axis equal
    grid on
    az = 135;
    el = 30;
    view(az,el)

    % Geometric reference frame
    x_g   = quiver3(0,0,0,1,0,0,color='r',LineWidth=1,MaxHeadSize=1);
    y_g   = quiver3(0,0,0,0,1,0,color='r',LineWidth=1,MaxHeadSize=1);
    z_g   = quiver3(0,0,0,0,0,1,color='r',LineWidth=1,MaxHeadSize=1);

    % Body reference frame
    %x_b   = quiver3(0,0,0,L_gb(1,1),L_gb(2,1),L_gb(3,1),color='g',LineWidth=1,MaxHeadSize=1);
    %y_b   = quiver3(0,0,0,L_gb(1,2),L_gb(2,2),L_gb(3,2),color='g',LineWidth=1,MaxHeadSize=1);
    %z_b   = quiver3(0,0,0,L_gb(1,3),L_gb(2,3),L_gb(3,3),color='g',LineWidth=1,MaxHeadSize=1);

    % LVLH reference frame
    x_LVLH   = quiver3(0,0,0,L_gLVLH(1,1),L_gLVLH(2,1),L_gLVLH(3,1),color='b',LineWidth=1,MaxHeadSize=1);
    y_LVLH   = quiver3(0,0,0,L_gLVLH(1,2),L_gLVLH(2,2),L_gLVLH(3,2),color='b',LineWidth=1,MaxHeadSize=1);
    z_LVLH   = quiver3(0,0,0,L_gLVLH(1,3),L_gLVLH(2,3),L_gLVLH(3,3),color='b',LineWidth=1,MaxHeadSize=1);

    % NED reference frame
    x_ned   = quiver3(0,0,0,L_gned(1,1),L_gned(2,1),L_gned(3,1),color='c',LineWidth=1,MaxHeadSize=1);
    y_ned   = quiver3(0,0,0,L_gned(1,2),L_gned(2,2),L_gned(3,2),color='c',LineWidth=1,MaxHeadSize=1);
    z_ned   = quiver3(0,0,0,L_gned(1,3),L_gned(2,3),L_gned(3,3),color='c',LineWidth=1,MaxHeadSize=1);


    % Labels
    text(1, 0, 0, 'X', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(0, 1, 0, 'Y', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    text(0, 0, 1, 'Z', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 7);
    %text(L_gb(1,1), L_gb(2,1), L_gb(3,1), 'X_b', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    %text(L_gb(1,2), L_gb(2,2), L_gb(3,2), 'Y_b', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    %text(L_gb(1,3), L_gb(2,3), L_gb(3,3), 'Z_b', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(L_gLVLH(1,1), L_gLVLH(2,1), L_gLVLH(3,1), 'X_{LVLH}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(L_gLVLH(1,2), L_gLVLH(2,2), L_gLVLH(3,2), 'Y_{LVLH}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(L_gLVLH(1,3), L_gLVLH(2,3), L_gLVLH(3,3), 'Z_{LVLH}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(L_gned(1,1), L_gned(2,1), L_gned(3,1), 'X_{NED}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(L_gned(1,2), L_gned(2,2), L_gned(3,2), 'Y_{NED}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);
    text(L_gned(1,3), L_gned(2,3), L_gned(3,3), 'Z_{NED}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7);

    % Velocities
    v_satg_vector   = quiver3(0,0,0,v_satg(1)/magnv_sat,v_satg(2)/magnv_sat,v_satg(3)/magnv_sat,color='k',LineWidth=2,MaxHeadSize=2); 
    v_satg_dir      = quiver3(0,0,0,v_satg(1)/magnv_sat*10,v_satg(2)/magnv_sat*10,v_satg(3)/magnv_sat*10,color='k',LineStyle='--',LineWidth=1,MaxHeadSize=0);
    v_windg_vector  = quiver3(0,0,0,v_windg(1)/magnv_sat,v_windg(2)/magnv_sat,v_windg(3)/magnv_sat,color='k',LineWidth=2,MaxHeadSize=2); 
    v_windg_dir     = quiver3(0,0,0,v_windg(1)/magnv_sat*1000,v_windg(2)/magnv_sat*1000,v_windg(3)/magnv_sat*1000,color='k',LineStyle='--',LineWidth=1,MaxHeadSize=0);
    v_relg_vector   = quiver3(0,0,0,v_relg(1)/magnv_sat,v_relg(2)/magnv_sat,v_relg(3)/magnv_sat,color='k',LineWidth=2,MaxHeadSize=2); 
    v_relg_dir      = quiver3(0,0,0,v_relg(1)/magnv_sat*10,v_relg(2)/magnv_sat*10,v_relg(3)/magnv_sat*10,color='k',LineStyle='--',LineWidth=1,MaxHeadSize=0); 

    axis equal
    xlim([-axlength,axlength])
    ylim([-axlength,axlength])
    zlim([-axlength,axlength])
    xlabel('x_{g}'); ylabel('y_{g}'); zlabel('z_{g}')
    









end

%------------- END OF CODE --------------
