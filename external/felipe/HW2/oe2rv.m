function [rv] = oe2rv(oe,deltat,mu)
    a     = oe(1);
    e     = oe(2);
    i     = oe(3);
    omega = oe(4);
    Omega = oe(5);
    M0    = oe(6);
    
    if deltat == 0
       M = M0;
    else
       M = wrapTo2Pi(M0 + deltat*sqrt(mu/(a^3)));
    end
    
    f  = @(E) E - e*sin(E) - M;
    E0 = M;
    E  = fzero(f,E0);
    nu = 2*atan2(sqrt(1 + e)*sin(E/2),sqrt(1 - e)*cos(E/2));
    
    ogoe = [a e i omega Omega nu];
    rv   = ogoe2rv(ogoe,mu);
end