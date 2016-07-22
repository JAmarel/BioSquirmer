function [v1, v2] = Rotate_Vector(v1, v2, theta)
%Uses the rotation matrix to rotate the input vectors by theta
v1New = v1*cos(theta) - v2*sin(theta);
v2New = v1*sin(theta) + v2*cos(theta);
v1 = v1New;
v2 = v2New;

end

