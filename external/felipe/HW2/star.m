function qstar = star(q)
    qstar      = zeros(4,1);
    qstar(1:3) = -q(1:3);
    qstar(4)   = q(4);
end