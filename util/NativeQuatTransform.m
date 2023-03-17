function v_B = NativeQuatTransform(q_A2B,v_A)
%NativeQuatTransform performs a coordinate transformation on vector, v_A, using
%quaternion q_A2B without converting q_A2B into a DCM.

v_B_prime = QuatProduct(q_A2B,QuatProduct([v_A; 0],QuatInv(q_A2B)));
v_B = v_B_prime(1:3);

end