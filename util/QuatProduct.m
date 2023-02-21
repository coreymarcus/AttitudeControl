function q_A2C = QuatProduct(q_B2C,q_A2B)
%QuatProduct performs the Shuster quaternion product with scalar last
%quaternions

% Extract components
s_A2B = q_A2B(4,1);
s_B2C = q_B2C(4,1);
v_A2B = q_A2B(1:3,1);
v_B2C = q_B2C(1:3,1);

% Multiplication
s_A2C = s_A2B*s_B2C - dot(v_A2B,v_B2C);
v_A2C = s_A2B*v_B2C + s_B2C*v_A2B - cross(v_B2C,v_A2B); % I think this order is actually correct
% v_A2C = s_A2B*v_B2C + s_B2C*v_A2B - cross(v_A2B,v_B2C);

% Output
q_A2C = [v_A2C; s_A2C];

end