function T = Quat2DCM(q)
%Quat2DCM converts a quaternion to passive rotation (transformation) matrix

% Convert to axis angle
[theta, n] = Quat2AxisAngle(q);

% Convert to transformation matrix
T = AxisAngle2DCM(theta,n);

end