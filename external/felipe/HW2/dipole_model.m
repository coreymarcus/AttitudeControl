function Bn = dipole_model(r,lat)
    B_0 = 3.12e-5;   % Tesla
    R_e = 6.3781e+6; % m
    Bn  = B_0*((R_e/norm(r))^3)*[cos(lat) 0 2*sin(lat)]';
end