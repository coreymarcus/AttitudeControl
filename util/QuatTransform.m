function v_B = QuatTransform(q_A2B,v_A)
%QuatTransform performs a coordinate transformation on vector, v_A, using
%quaternion q_A2B.
v_B = Quat2DCM(q_A2B)*v_A;

end