function [ pathOut ] = ADBSatFcn_eulerAngles( modpath, respath, param_eq, yaw_deg, pitch_deg, roll_deg, phi_deg, inc_deg,RAAN,omega,theta, magn_vsat, v_corot, flag_shadow, flag_solar,flag_wind,flag_corot, env, del, verb)
%ADBSATFCN Creates a .mat file(s) in "/inou/results" with the following fields:
%
% Inputs:
%       modname : name of the input model file (.mat) in /inou/models
%       param_eq : GSI & SRP parameter structure inputs (dependent on GSI model)
%           gsi_model   : GSI Model (string)
%           alpha       : Energy accommodation coefficient (Cook/Sentman/Maxwell)
%           alphaN      : Normal thermal accommodation coefficient (CLL)
%           sigmaN      : Normal momentum accommodation coefficent (S&C/Storch)
%           sigmaT      : Tangential momentum accommodation coefficient (S&C/Storch/CLL)
%           sol_cR      : Solar specular reflectivity coefficient
%           sol_cD      : Solar diffuse reflectivity coefficient
%       yaw_deg   : vector of yaw angles [deg]
%       pitch_deg : vector of pitch angles [deg]
%       roll_deg  : vector of roll angles [deg]
%       phi_deg   : angle between flight direction and position vector [deg]
%       inc_deg   : inclination of orbit [deg]
%       RAAN      : Right-ascension of Ascending Node [rad]
%       omega     : Argument of Perigee [rad]
%       theta     : True Anomaly [rad]
%       magn_vsat : Magnitude of satellite track velocity [m/s]
%       v_corot   : Induced velocity through corrotation of atmosphere [m/s]
%       shadow    : flag for shadown analysis
%       solar     : flag for solar coefficient analysis
%       wind      : flag for consideration of thermospheric wind effects
%       corot     : flag for consideration of the co-rotation of the atmosphere
%       env       : vector of input environmental parameters
%       del       : flag for individual file cleanup (on mergeAEDB)
%       verb      : flag for visual output
%
% Outputs:
%       fileout  : .mat file(s) with the following items in a structure variable
%           yaw         : Yaw angle
%           pitch       : Pitch angle
%           roll        : Roll angle
%           Aref        : Reference area = Area total/2
%           Lref        : Reference longitude = Xmax-Xmin
%           AreaProj    : Projected area
%           aero
%               Cf_w     : Total aerodynamic force coeffients in wind axes [x,y,z] (3x1)
%               Cf_f     : Total aerodynamic force coeffients in flight axes [x,y,z] (3x1)
%               Cm_B     : Total aerodynamic moment coeffients in body axes [x,y,z](3x1)
%           solar
%               C_s      : Total solar force coefficients in
%               C_sB     : Total solar moment coefficients in
%           param_eq : Structure containing input GSI parameters
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

% Check GSIM Inputs
if strcmpi(param_eq.gsi_model,'CLL')
    if any(~isfield(param_eq,{'alphaN','sigmaT'}))
        error([param_eq.gsi_model,' model requires alphaN & sigmaT inputs'])
    end
elseif any(strcmpi(param_eq.gsi_model,{'schaaf','storchHyp'}))
    if any(~isfield(param_eq,{'sigmaN','sigmaT'}))
        error([param_eq.gsi_model,' model requires sigmaN and sigmaT inputs'])
    end
elseif any(strcmpi(param_eq.gsi_model,{'cook','sentman','maxwell','DRIA'}))
    if ~isfield(param_eq,'alpha')
        error([param_eq.gsi_model,' model requires alpha input'])
    end
elseif strcmpi(param_eq.gsi_model,'newton')
    % Nothing
else
    error('GSI Model input string not recognised');
end

% Check SRP Model Inputs
if flag_solar == 1 && any(~isfield(param_eq,{'sol_cR','sol_cD'}))
    error('SRP interaction model requires sol_cR and sol_cD inputs')
end

if ~isfield(param_eq,'Tw')
    warning('Wall temperature not defined in input. Default value of 300K used')
    param_eq.Tw = 300; % Surface temperature [K]
end


if numel(env) > 1
    % Environment Calculations
    param_eq = environment(param_eq, env(1),env(2),env(3),env(4),env(5),env(6),env(7),env(8),env(9:15),env(16));
end

%param_eq.Vw = sqrt(pi.*Rmean.*param_eq.Tw/2); % Average velocity of the reflected diffuse molecules

% Calculating the horizontal wind velocity
if flag_wind
    v_wind = atmoshwm(env(2),env(3),env(1),'day',env(5),'seconds',env(6),'model','total', 'version', '14'); % Horizontal wind in m/s
    v_windmzu = [transpose(v_wind);0]; % Wind velocity in Meridian-Zonal-Up Frame
else
    v_windmzu = [0;0;0]; % Wind velocity in Meridian-Zonal-Up Frame
end

% Setting the velocity induced by co-rotation
if flag_corot
    v_coroteic = v_corot;
else
    v_coroteic = [0;0;0];
end

% Calculate Interactions
pathOut = calc_coeff_eulerAngles(   modpath, respath, ...
                                    deg2rad(yaw_deg), deg2rad(pitch_deg), deg2rad(roll_deg), deg2rad(phi_deg), deg2rad(inc_deg), ...
                                    RAAN,omega,theta,v_windmzu,magn_vsat,v_coroteic,param_eq, flag_shadow, flag_solar, del, verb);

%------------- END OF CODE --------------
