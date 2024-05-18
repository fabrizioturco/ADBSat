function [ADBout] = SinglePoint_inOrbit_Analysis(alt,inc,lat,lon,year,dayofyear,UTseconds,gravParam,phi,magn_vsat,v_corot,kep,flags)

%% SETTINGS
modName = 'BEESAT4_lowFidelity_4';

% Path to model file
ADBSat_path = ADBSat_dynpath;
modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);
modOut = fullfile(ADBSat_path,'inou','models');
resOut = fullfile(ADBSat_path,'inou','results');

%% INPUT
%Input conditions
mass = 1.030;
f107Average = 140; % 65, 140, 250
f107Daily = 140; % 65, 140, 250
magneticIndex = ones(1,7)*15; % 0, 15, 45

% Initialize flags
verb = flags(1);
del = flags(2);
shadow = flags(3);
AnO = flags(4);
solar = flags(5);
wind = flags(6);
corot = flags(7);

% Initialize needed keplerian elements
RAAN = kep(4);
omega = kep(5);
theta = kep(6);
  
env = [alt*1e3, lat, lon, year, dayofyear, UTseconds, f107Average, f107Daily, magneticIndex, AnO, gravParam]; % Environment variables

% Attitude
yaw     = 0;  % Yaw angle [deg]
pitch   = 0;  % Pitch angle [deg]
roll    = 0; % Roll angle [deg]


% Model parameters
inparam.gsi_model = 'sentman';
% inparam.alpha = 0.85; % Accommodation (altitude dependent)
% inparam.sigmaN = 0.85;
% inparam.sigmaT = 0;

inparam.Tw = 300; % Wall Temperature [K]

inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

if AnO
    Oflag = 'Oxygen';
else
    Oflag = 'NoOxygen';
end

% SESAM: Semiempirical Model for Satellite Energy-Accommodation Coefficients
[T, rho] = atmosnrlmsise00(alt*1e3, lat, lon, year, dayofyear, UTseconds, f107Average, f107Daily, magneticIndex, Oflag);
inparam.alpha = 7.5E-17*rho(2)*T(2) / (1+7.5E-17*rho(2)*T(2)); % SESAM for accommodation coefficient (altitude dependent)

%% ANALYSIS
tic
% Import model
[modOut] = ADBSatImport(modIn, modOut, 0);

% Calculate
[ADBout] = ADBSatFcn_eulerAngles(modOut, resOut, inparam, yaw, pitch, roll, phi, inc, RAAN, omega, theta, magn_vsat, v_corot, shadow, solar, wind, corot, env, del, verb);
result = load(ADBout);

time = toc;
%% RESULTS
% Plot
if verb && ~del
    az = 135;
    el = 30;
    plot_surfq_LVLH(ADBout, modOut, yaw(1)*pi/180, pitch(1)*pi/180, roll(1)*pi/180, 'cp',1,az,el);
end

fprintf('\n')
% fprintf('C_D \t= %.4g\n', -result.Cf_w(1));
fprintf('C_D \t= %.4g\n', -result.Cf_w(1)*result.Aref/result.AreaProj);
fprintf('beta \t= %.4g\n', (result.Aref*(-result.Cf_w(1)))/mass);
fprintf('alpha \t= %.4g\n', inparam.alpha);
fprintf('A_ref \t= %.4g m^2\n', result.Aref);
fprintf('A_proj \t= %.4g m^2\n', result.AreaProj);
fprintf('L_ref \t= %.4g m\n', result.Lref);
fprintf('rho \t= %.4g kg/m^3\n', rho(6));
fprintf('T \t\t= %.4g K\n', T(:));
fprintf('\nComp time \t= %.4g s\n', time);
%% 
fprintf('\n\n')
end
%------------ END CODE -----------%

