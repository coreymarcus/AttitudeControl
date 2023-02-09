function qinv = QuatInv(q)
%QuatInv provides the inverse of a scalar last quaternion. Quaternion is a
%column vector
qinv = [-1*q(1:3); q(4)];
end