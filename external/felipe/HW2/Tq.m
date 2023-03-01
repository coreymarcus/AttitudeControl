function T = Tq(q)
    qv = q(1:3);
    qs = q(4);

    qcross = [0 -qv(3) qv(2); qv(3) 0 -qv(1); -qv(2) qv(1) 0];

    T = eye(3) - 2*qs*qcross + 2*qcross*qcross;
end

