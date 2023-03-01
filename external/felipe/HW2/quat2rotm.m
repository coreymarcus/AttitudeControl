function T = quat2rotm(q)

    qvec = q(1:3);
    q4   = q(4);
   
    q1 = qvec(1);
    q2 = qvec(2);
    q3 = qvec(3);
    T  = [q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4); ...
          2*(q2*q1-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4); ...
          2*(q3*q1+q2*q4) 2*(q3*q2-q1*q4) -q1^2-q2^2+q3^2+q4^2];
end