function [theta,n] = DCM2AxisAngle(T)
%DCM2AxisAngle

% Find eigenvalues of T
[V, D] = eig(T);
lambdas = diag(D);
targ = lambdas == 1;

% Make sure we've found a correct eigenvalue
if(sum(targ) ~= 1)
    error("No unity eigenvalue found.")
end

% Eigenvector is the rotation axis
n = V(:,targ);
assert(norm(n) == 1);

% Angle is function of the trace
theta = acos(0.5*trace(T) - 0.5);

% Check to verify we have the right angle
T2 = AxisAngle2DCM(theta,n);
if(norm(T2 - T,"fro") > 0)
    theta = -1*theta;
end

end