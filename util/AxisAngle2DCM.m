function T = AxisAngle2DCM(theta,n)
%AxisAngle2DCM is self explanitory

nx = CrossProductMat(n);
T = eye(3) + sin(theta)*nx + (1 - cos(theta))*nx*nx;

end