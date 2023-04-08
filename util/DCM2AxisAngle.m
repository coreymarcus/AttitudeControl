function [theta,n] = DCM2AxisAngle(T)
%DCM2AxisAngle

% Find eigenvalues of T
[V, D] = eig(T);
lambdas = diag(D);
targ = find(abs(real(lambdas) - 1) < 1E-8); % This is inefficient but doing it this way for simulink compatibility

% Check for no rotation
if(length(targ) == 3)
    theta = 0;
    n = [1 0 0]';
    return;
end

% Make sure we've found a correct eigenvalue
if(length(targ) ~= 1)
    n = [1 0 0]';
    theta = 0;
    error("No unity eigenvalue found.")
    return;
end

% Eigenvector is the rotation axis
n = real(V(1:3,targ(1)));
assert((norm(n) - 1) < 1E-8);
if(norm(imag(V(1:3,targ(1)))) > 0)
    error("Imaginary Rotation Axis")
end

% Angle is function of the trace
theta = acos(0.5*trace(T) - 0.5);

% Check to verify we have the right angle
T2 = AxisAngle2DCM(theta,n);
if(norm(T2 - T,"fro") > norm(T2' - T,"fro"))
    theta = -1*theta;
end

end