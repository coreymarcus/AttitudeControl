function [theta, n] = Quat2AxisAngle(q)
%Quat2AxisAngle converts quaternion to axis angle representation.
%Quaternions are scalar last and right handed

% if(abs(q(4)) > 1)
%     x = 5;
% end

theta = 2*acos(q(4));

if(theta > 0)
    n = q(1:3)/sin(theta/2);
else
    % Need to handle zero rotation case
    n = [1 1 1]/sqrt(3);
end

end