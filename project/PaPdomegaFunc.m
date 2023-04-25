function B = PaPdomegaFunc(dt,w_hat1,w_hat2,w_hat3)
%PaPdomegaFunc
%    B = PaPdomegaFunc(DT,W_HAT1,W_HAT2,W_HAT3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    24-Apr-2023 20:57:34

et1 = integral(@(tau)(w_hat2.^2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)+(w_hat3.^2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)+1.0,0.0,dt);
t2 = et1;
t3 = integral(@(tau)w_hat3.*sin(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)).*(dt-tau).*1.0./sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)-(w_hat1.*w_hat2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2),0.0,dt);
t4 = integral(@(tau)-w_hat2.*sin(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)).*(dt-tau).*1.0./sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)-(w_hat1.*w_hat3.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2),0.0,dt);
t5 = integral(@(tau)-w_hat3.*sin(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)).*(dt-tau).*1.0./sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)-(w_hat1.*w_hat2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2),0.0,dt);
et2 = integral(@(tau)(w_hat1.^2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)+(w_hat3.^2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)+1.0,0.0,dt);
t6 = et2;
t7 = integral(@(tau)w_hat1.*sin(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)).*(dt-tau).*1.0./sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)-(w_hat2.*w_hat3.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2),0.0,dt);
t8 = integral(@(tau)w_hat2.*sin(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)).*(dt-tau).*1.0./sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)-(w_hat1.*w_hat3.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2),0.0,dt);
t9 = integral(@(tau)-w_hat1.*sin(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)).*(dt-tau).*1.0./sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)-(w_hat2.*w_hat3.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2),0.0,dt);
et3 = integral(@(tau)(w_hat1.^2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)+(w_hat2.^2.*(cos(sqrt(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2))-1.0).*(dt-tau).^2)./(w_hat1.^2.*(dt-tau).^2+w_hat2.^2.*(dt-tau).^2+w_hat3.^2.*(dt-tau).^2)+1.0,0.0,dt);
t10 = et3;
B = reshape([t2,t5,t8,t3,t6,t9,t4,t7,t10],[3,3]);
