function [x] = Keplerian2ECI(kep,mu)

%% Determining r and v from the orbital elements %%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         
% This function converts keplerian elements to a position vector expressed
% in the Earth Centered Inertial reference system (ECI) trough coordinate
% transformation.
%
%                                                                         
%  Author: Laurin Mächtig
%
%
% Date: 01.05.2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Divide input vector in singular keplerian elements %%%%%%%%%%%%%%%%%%%%%

a  = kep(1);
e  = kep(2);
i  = kep(3);
RAAN  = kep(4);
omega  = kep(5);
theta = kep(6);

%% Express state vector in the Perifocal reference frame %%%%%%%%%%%%%%%%%%

% Perifocal reference frame (origin: earth; x: points to perigee, z: parallel to h)
p = a*(1-e^2);
r_perifocal = p/(1+e*cos(theta))*[cos(theta);sin(theta);0];
v_perifocal = sqrt(mu/p)*[-sin(theta);e+cos(theta);0];

%% Conversion to ECI (Earth centered intertial reference frame) %%%%%%%%%%%

rotate_xp2node = angle2dcm(-omega,0,0);                                     % rotate around z-axis until x_perifocal is aligned with the vector pointing to the ascending node
rotate_inc   = angle2dcm(0,0,-i);                                           % rotate around new x-axis in the equatorial plane until the z-axis aligns with the Intertial reference frame z-axis
rotate_xnew2xeic = angle2dcm(-RAAN,0,0);                                    % rotate around new z-axis until x-axis is aligned with the Inertial reference frame x-axis
L_eicPerifocal = rotate_xnew2xeic * rotate_inc * rotate_xp2node;

r_eic = L_eicPerifocal * r_perifocal;
v_eic = L_eicPerifocal * v_perifocal;

%% St Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [r_eic; v_eic];

end
