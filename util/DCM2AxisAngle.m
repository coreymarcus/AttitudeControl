function [theta,n] = DCM2AxisAngle(T)
%DCM2AxisAngle

% Find eigenvalues of T
[V, D] = eig(T);
lambdas = diag(D);
targ = lambdas == 1;

% Check for no rotation
if(sum(targ) == 3)
    theta = 0;
    n = [1 0 0]';
    return;
end

% Make sure we've found a correct eigenvalue
if(sum(targ) ~= 1)
    n = [1 0 0]';
    theta = 0;
    error("No unity eigenvalue found.")
    return;
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