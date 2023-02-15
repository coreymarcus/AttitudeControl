function q3 = schuster_mult(q1,q2)
    qv1 = q1(1:3);
    qs1 = q1(4);

    qv2 = q2(1:3);
    qs2 = q2(4);

    qv3 = qs2*qv1 + qs1*qv2 - cross(qv1,qv2);
    qs3 = qs1*qs2 - dot(qv1,qv2);

    q3 = [qv3;qs3];
end