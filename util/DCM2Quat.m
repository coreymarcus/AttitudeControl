function q = DCM2Quat(T)
%DCM2Quat converts

[theta,n] = DCM2AxisAngle(T);
q = AxisAngle2Quat(theta,n);

end