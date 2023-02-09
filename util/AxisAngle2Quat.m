function q = AxisAngle2Quat(theta,n)
%AxisAngle2Quat provides the conversion
q = [n*sin(theta/2);
    cos(theta/2)];
end