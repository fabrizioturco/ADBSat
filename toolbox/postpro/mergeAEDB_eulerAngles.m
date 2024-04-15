function fileOut = mergeAEDB_eulerAngles(pathName, fiName, del)
% Merges coefficients computed for different aoa/aos into a single structure
% Assumes that a complete grid of aoa and aos values is available in the
% input folder
%
% Inputs:
%       pathName        : Path of the folder containing the multiple .mat files
%       fiName          : Name of the model file
%
% Output: (produces one output file containing a structure with the following fields)
%        yaw         : Yaw angles (rad)
%        pitch       : Pitch angles (rad)
%        roll        : Roll angles (rad)
%        Aref     : Reference area
%        Lref      : Refrence length
%        AreaProj    : Projected areas
%        aero
%            Cf_w    : Aerodynamic force coeffients in wind axes [x,y,z] (yaw x pitch x roll)
%            Cf_f    : Aerodynamic force coeffients in flight axes [x,y,z] (yaw x pitch x roll)
%            Cm_B    : Aerodynamic moment coeffients in body axes [x,y,z] (yaw x pitch x roll)
%        solar
%            C_s     : Total solar force coefficients in axes [x,y,z] (yaw x pitch x roll)
%            C_sB    : Total solar moment coefficients in axes [x,y,z] (yaw x pitch x roll)
%        param_eq    : Structure containing input GSI parameters
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
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

%Filestructure
fis = dir(pathName);
for ii = 1:length(fis)
    ismod(ii) = strncmp(fis(ii).name,fiName,length(fiName));
end
nFiles = sum(ismod);

% Preallocate
v_yaw   = zeros(1,nFiles);
v_pitch = zeros(1,nFiles);
v_roll  = zeros(1,nFiles);
v_Cf_w  = zeros(3,nFiles);
v_Cf_f  = zeros(3,nFiles);
v_Cm_B  = zeros(3,nFiles);
v_AreaProj = zeros(1,nFiles);
v_Cf_s  = zeros(3,nFiles);
v_Cm_S = zeros(3,nFiles);
check_input = 1;

% Loop through each matching file and load in parameters
count = 0;
for ii = 1:length(fis)
    if ~ismod(ii)
        continue
    else
        load(fullfile(pathName,fis(ii).name));
        count = count+1;
        
        Aref_old = Aref;
        Lref_old = Lref;
        param_eq_old = param_eq;        
        
        v_yaw(count)        = yaw;
        v_pitch(count)      = pitch;
        v_roll(count)       = roll;
        v_Cf_w(:,count)     = Cf_w;
        v_Cf_f(:,count)     = Cf_f;
        v_Cm_B(:,count)     = Cm_B;
        v_AreaProj(count)   = AreaProj;
        
        if exist('Cf_s','var')
            v_Cf_s(:,count)  = Cf_s;
            v_Cm_S(:,count)  = Cm_S;
        else
            v_Cf_s(:,count)  = NaN;
            v_Cm_S(:,count)  = NaN;
        end
        
        if ~isequal(param_eq, param_eq_old) || ~isequal(Aref, Aref_old) || ~isequal(Lref, Lref_old)
            warning('Change in input model or GSI paramaters detected...')
            check_input = 0;
        end
        
        clear yaw pitch roll Cf_w Cf_f Cm_B AreaProj Cf_s Cm_S
    end
end
[~,index] = sortrows([v_yaw', v_pitch', v_roll']);
[yaw_u, ~] = unique(v_yaw); % Unique Yaw
[pitch_u, ~] = unique(v_pitch); % Unique Pitch
[roll_u, ~] = unique(v_roll); % Unique Roll

n_yaw   = numel(yaw_u);
n_pitch = numel(pitch_u);
n_roll  = numel(roll_u);

[aedb.yaw, aedb.pitch, aedb.roll] = meshgrid(yaw_u, pitch_u, roll_u);

aedb.aero.Cf_wX = reshape(v_Cf_w(1,index), [n_yaw n_pitch n_roll])';
aedb.aero.Cf_wY = reshape(v_Cf_w(2,index), [n_yaw n_pitch n_roll])';
aedb.aero.Cf_wZ = reshape(v_Cf_w(3,index), [n_yaw n_pitch n_roll])';

aedb.aero.Cf_fX = reshape(v_Cf_f(1,index), [n_yaw n_pitch n_roll])';
aedb.aero.Cf_fY = reshape(v_Cf_f(2,index), [n_yaw n_pitch n_roll])';
aedb.aero.Cf_fZ = reshape(v_Cf_f(3,index), [n_yaw n_pitch n_roll])';

aedb.aero.Cm_BX = reshape(v_Cm_B(1,index), [n_yaw n_pitch n_roll])';
aedb.aero.Cm_BY = reshape(v_Cm_B(2,index), [n_yaw n_pitch n_roll])';
aedb.aero.Cm_BZ = reshape(v_Cm_B(3,index), [n_yaw n_pitch n_roll])';

aedb.AreaProj = reshape(v_AreaProj(1,index), [n_yaw n_pitch n_roll])';

aedb.solar.Cf_sX = reshape(v_Cf_s(1,index), [n_yaw n_pitch n_roll])';
aedb.solar.Cf_sY = reshape(v_Cf_s(2,index), [n_yaw n_pitch n_roll])';
aedb.solar.Cf_sZ = reshape(v_Cf_s(3,index), [n_yaw n_pitch n_roll])';

aedb.solar.Cm_sBX = reshape(v_Cm_S(1,index), [n_yaw n_pitch n_roll])';
aedb.solar.Cm_sBY = reshape(v_Cm_S(2,index), [n_yaw n_pitch n_roll])';
aedb.solar.Cm_sBZ = reshape(v_Cm_S(3,index), [n_yaw n_pitch n_roll])';

if check_input
    aedb.Aref = Aref;
    aedb.Lref = Lref;
    aedb.param_eq = param_eq;
end

fileOut = fullfile([pathName,'.mat']);
save(fileOut, 'aedb');

% Delete individual files and folder
if del
    rmdir(pathName,'s')
end

%------------- END OF CODE --------------
