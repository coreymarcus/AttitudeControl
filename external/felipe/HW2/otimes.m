function w = otimes(qbar,q)
    qbarv = qbar(1:3);
    qbars = qbar(4);
    qv    = q(1:3);
    qs    = q(4);
    wv    = qs*qbarv + qbars*qv - cross(qbarv,qv);
    ws    = qbars*qs - dot(qbarv,qv);
    w     = [wv;ws];
end