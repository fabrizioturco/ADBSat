clear
close all

%% SETTINGS
% modName = 'BEESAT4_highFidelity'
% modName = 'BEESAT4_highFidelity_noAntennae';
% modName = 'BEESAT4_highFidelity_flatPlate';
% modName = 'BEESAT4_lowFidelity_3';
% modName = 'BEESAT4_lowFidelity_3_noAntennae';
% modName = 'BEESAT4_lowFidelity_2';
% modName = 'BEESAT4_lowFidelity_2_noAntennae';
modName = 'BEESAT4_lowFidelity_4';
% modName = 'BEESAT4_lowFidelity_noAntennae';



% Path to model file
ADBSat_path = ADBSat_dynpath;
modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);
modOut = fullfile(ADBSat_path,'inou','models');
resOut = fullfile(ADBSat_path,'inou','results');

%% INPUT
%Input conditions
mass = 1.030;
alt = 400; %km
inc = 97.03; %deg
lat = 0;
lon =  0;
year = 2023;
dayofyear = 185; % 04.07.2023
UTseconds = 12*3600;%16*3600 + 29*60 + 32;
f107Average = 140; % 65, 140, 250
f107Daily = 140; % 65, 140, 250
magneticIndex = ones(1,7)*15; % 0, 15, 45
AnO =  1;   
gravParam = 3.986*10^14; % m^3/s^2
env = [alt*1e3, lat, lon, year, dayofyear, UTseconds, f107Average, f107Daily, magneticIndex, AnO, gravParam]; % Environment variables

% Attitude
yaw     = 0;  % Yaw angle [deg]
pitch   = 0;  % Pitch angle [deg]
roll    = 0; % Roll angle [deg]

phi = 60; % Angle between flight direction and position vector [deg]


% Model parameters
shadow = 1;
inparam.gsi_model = 'sentman';
% inparam.alpha = 0.85; % Accommodation (altitude dependent)
% inparam.sigmaN = 0.85;
% inparam.sigmaT = 0;

inparam.Tw = 300; % Wall Temperature [K]

solar = 0;
inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

wind = 1; % Flag for consideration of thermospheric wind effects

verb = 1;
del = 0;

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
[ADBout] = ADBSatFcn_eulerAngles(modOut, resOut, inparam, yaw, pitch, roll, phi, inc, shadow, solar,wind, env, del, verb);
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
%------------ END CODE -----------%