function T = Quat2DCM(q)
%Quat2DCM converts a quaternion to passive rotation (transformation) matrix

% Convert to axis angle
[theta, n] = Quat2AxisAngle(q);

% Convert to transformation matrix
nx = CrossProductMat(n);
T = eye(3) + sin(theta)*nx + (1 - cos(theta))*nx*nx;

end