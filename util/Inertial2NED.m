function T_inertial2NED = Inertial2NED(r_inertial)
%Inertial2NED generates the transformation between inertial and NED
%coordinates as a function of the position in inertial (earth centered)
%coordinates

% Inertial z axis
iz = [0 0 1]';

% componenets
nz = -r_inertial/norm(r_inertial);
ny = cross(nz,iz)/norm(cross(nz,iz));
nx = cross(ny,nz);

% Build matrix
T_inertial2NED = [nx';
    ny';
    nz'];

end