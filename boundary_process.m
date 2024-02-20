function [A,b] = boundary_process(A,b)
A(1,:) = 0;
A(end,:) = 0;
A(1,1) = 1;
A(end,end) = 1;
b(1) = 0;
b(end) = cos(1);
end