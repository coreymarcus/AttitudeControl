function Tn_i = Tn_i(r)
    ni_z =  - r/norm(r);
    ni_y = cross(ni_z,[0 0 1]')/norm(cross(ni_z,[0 0 1]'));
    ni_x = cross(ni_y,ni_z);

    Tn_i = [ni_x'; ni_y'; ni_z'];
end