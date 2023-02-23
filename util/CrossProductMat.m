function out = CrossProductMat(vec)
%CROSSPRODUCTMAT creates a cross product equivalent matrix
arguments
    vec (3,1)
end 
out = [0 -vec(3) vec(2);
    vec(3) 0 -vec(1);
    -vec(2) vec(1) 0];

end

