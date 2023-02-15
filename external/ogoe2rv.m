function [rv] = ogoe2rv(oe,mu)
    a     = oe(1);
    e     = oe(2);
    i     = oe(3);
    omega = oe(4);
    Omega = oe(5);
    nu    = oe(6);
    
    p = a*(1-e^2);
    r = p/(1+e*cos(nu));
    
    rpqw = [r*cos(nu);r*sin(nu);0];
    vpqw = [-sqrt(mu/p)*sin(nu);sqrt(mu/p)*(e+cos(nu));0];
    
    r = rotate(-Omega,3)*rotate(-i,1)*rotate(-omega,3)*rpqw;
    v = rotate(-Omega,3)*rotate(-i,1)*rotate(-omega,3)*vpqw;
    
    rv = [r;v];
end