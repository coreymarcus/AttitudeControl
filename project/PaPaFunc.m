function A = PaPaFunc(dt,w_hat1,w_hat2,w_hat3)
%PaPaFunc
%    A = PaPaFunc(DT,W_HAT1,W_HAT2,W_HAT3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    27-Jun-2023 14:03:21

t2 = dt.^2;
t3 = w_hat1.^2;
t4 = w_hat2.^2;
t5 = w_hat3.^2;
t6 = t2.*t3;
t7 = t2.*t4;
t8 = t2.*t5;
t9 = t6+t7+t8;
t10 = 1.0./t9;
t11 = sqrt(t9);
t12 = cos(t11);
t13 = sin(t11);
t14 = 1.0./t11;
t15 = t12-1.0;
t16 = dt.*t13.*t14.*w_hat1;
t17 = dt.*t13.*t14.*w_hat2;
t18 = dt.*t13.*t14.*w_hat3;
t19 = t2.*t10.*t15.*w_hat1.*w_hat2;
t20 = t2.*t10.*t15.*w_hat1.*w_hat3;
t21 = t2.*t10.*t15.*w_hat2.*w_hat3;
t22 = t6.*t10.*t15;
t23 = t7.*t10.*t15;
t24 = t8.*t10.*t15;
t25 = -t19;
t26 = -t20;
t27 = -t21;
A = reshape([t23+t24+1.0,-t18+t25,t17+t26,t18+t25,t22+t24+1.0,-t16+t27,-t17+t26,t16+t27,t22+t23+1.0],[3,3]);
