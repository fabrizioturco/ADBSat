function [yawArray,pitchArray] = equidistributedEulerAngles(nPoints)
%EQUIDISTRIBUTEDEULERANGLES Creates arrays of yaw and pitch angles, so that flow direction relative to satellite is
%spherically evenly distributed.
% 
% Inputs:
%       nPoints : number off sampling points
%
% Outputs:
%       yawArray    : Yaw angles [deg]
%       pitchArray  : Pitch angles [deg]
%
% Author: Fabrizio Turco
% University of Stuttgart, Institute of Space Systems
% February 2024
% Reference: https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
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

if nargin == 0
    nPoints = 100;
end

r = 1;
a = 4 *pi * r^2/nPoints;
d = sqrt(a);
Mtheta= round(pi/d);
dtheta = pi/Mtheta;
dphi = a/dtheta;

Ncount = 0;
for m = 0:(Mtheta-1)
    theta = pi*(m + 0.5)/Mtheta;
    Mphi = round(2*pi*sin(theta)/dphi);
    for n = 0:(Mphi-1)
        thetaArray(Ncount+1) = theta;
        phiArray(Ncount+1) = 2*pi*n/Mphi;
        X(Ncount+1) = r * sin(thetaArray(Ncount+1)) * cos(phiArray(Ncount+1));
        Y(Ncount+1) = r * sin(thetaArray(Ncount+1)) * sin(phiArray(Ncount+1));
        Z(Ncount+1) = r * cos(thetaArray(Ncount+1));
        Ncount = Ncount + 1;
    end
end

pitchArray  = rad2deg(atan2(X,Z))';
yawArray    = rad2deg(asin(Y))';

end
%------------- END CODE --------------