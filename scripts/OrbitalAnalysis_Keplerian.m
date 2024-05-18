clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script initializes the satellite position through use of keplerian
% elements. These elements are then transformed into coordinates of the 
% earth centered inertial reference frame (ECI).
% Finally, the script calls ADBSat at mutliple positions of the chosen
% orbit to perform multiple analyses.
%
%
% Author: Laurin Mächtig
%
%
% Date: 01.05.2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Output path

currentDir = pwd;
inouDir = [currentDir, '\inou']; 
OutDir = [inouDir, '\Results_OrbitalAnalysis'];

if exist(OutDir, 'dir')
    rmdir(OutDir, 's');  % Remove old data
end
mkdir(OutDir);

%% Set initial keplerian elements

% Set constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 3.986e14;                                                              % Gravitational constant [m^3/s^2]
re  = 6378137;                                                              % Equatorial radius [m]
we  = 7.2921151 * 10^(-5);                                                  % Earth's rotation rate [rad/s]
we_vec = [0;0;we];

% Initialize keplerian elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a0       = re + 200 * 10^3;                                                 % Semi-major axis in m
e0       = 0.0001;                                                          % Eccentricity
i0       = deg2rad(98);                                                     % Inclination [rad]
RAAN0    = deg2rad(10);                                                     % Right Ascension of Ascending Node [rad]
omega0   = deg2rad(10);                                                     % Argument of perigee [rad]
theta0   = deg2rad(10);                                                     % True anomaly [rad]

kep0 = [a0; e0; i0; RAAN0; omega0; theta0];

% Loop over orbit positions
nPositions = 10;
for Pos = 1:nPositions
    theta = theta0 + 2*pi/nPositions * Pos;
    kep = [a0; e0; i0; RAAN0; omega0; theta];

%% Convert keplerian elements to state vector in ECI %%%%%%%%%%%%%%%%%%%%%%

    [x0] = Keplerian2ECI(kep,mu);
    r0 = x0(1:3);
    v0 = x0(4:6);

%% Calculate initialization variables for ADBSat %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    magn_vsat = sqrt(-mu/a0+2*mu/sqrt(dot(r0,r0)));                         % Satellite track velocity    

    v_corot = cross(we_vec,r0);                                            % Induced velocity through corrotation of the atmosphere

    inc = rad2deg(i0);                                                      % Inclination [deg]

    % Time specification
    year = 2023;
    month = 7;
    dayofmonth = 4;
    dayofyear = 185;                                                         % 04.07.2023
    hourofday = 12;
    minofhour = 0;
    secofmin = 0;
    UTseconds = hourofday*3600 + minofhour * 60 + secofmin;
    
    % Calculate Latitude and longitude from ECI coordinates %%%%%%%%%%%%%%%%%%%

    utc = [year dayofmonth month hourofday minofhour secofmin];             % Coordinated Universal Time
    lla = eci2lla(transpose(r0),utc);                                       % Transform position from ECI to LLA
    lat = deg2rad(lla(1));                                                  % Latitude [rad]
    lon = deg2rad(lla(2));                                                  % [Longitude [rad]
    alt = lla(3)/1000;                                                      % Altitude [km]

    
%% Calculate the angle between position and velocity vector %%%%%%%%%%%%%%%
    
    phi = rad2deg(acos(dot(r0,v0)/(norm(r0)*norm(v0))));                    % Angle between flight direction and position vector [deg]
    
%% Call ADBSat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set flags for analysis
    flag_verb = 0;
    flag_del = 0;
    flag_shadow = 1;
    flag_AnO = 1;
    flag_solar = 0;
    flag_wind = 1;
    flag_corot = 1;
    flags = [flag_verb, flag_del, flag_shadow, flag_AnO, flag_solar, flag_wind, flag_corot];
    
    [ADBout] = SinglePoint_inOrbit_Analysis(alt,inc,lat,lon,year,dayofyear,UTseconds,mu,phi,magn_vsat,v_corot,kep,flags);

%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    outName = "results_Pos" + Pos + ".mat";
    copyfile(ADBout, fullfile(OutDir,outName));
    delete(ADBout);

end

%% Calculate orbit averaged results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cfw_vec = zeros(nPositions,3);
C_D_vec = zeros(nPositions,1);
Aproj_vec = zeros(nPositions,1);
for Pos = 1:nPositions
    outName = "results_Pos" + Pos + ".mat";
    resStruct = load(fullfile(OutDir,outName));
    Cfw_vec(Pos,:) = resStruct.Cf_w;
    C_D_vec(Pos,1) = resStruct.C_D;
    Aproj_vec(Pos,1) = resStruct.AreaProj;
end

Cfw_averaged = [sum(Cfw_vec(:,1),"all"),sum(Cfw_vec(:,2),"all"),sum(Cfw_vec(:,3),"all")]/nPositions;
C_D_averaged = sum(C_D_vec,"all")/nPositions;
Aproj_averaged = sum(Aproj_vec,"all")/nPositions;

fprintf("Orbit-averaged values:\n")
fprintf('C_D \t= %.4g\n', C_D_averaged);
fprintf('A_proj \t= %.4g\n', Aproj_averaged);


%% 
fprintf('\n\n')
%------------ END CODE -----------%