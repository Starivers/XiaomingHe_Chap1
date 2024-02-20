function [Pb,Tb] = generate_PbTb(number_of_element,basis_type)

N = number_of_element;
switch basis_type
    case 101   % linear element, as well as mesh.
        Pb = 0:1/N:1; % default interval is [0,1], modify this part to change interval.
        Tb = [1:N; 2:N+1];
    case 102   % quadratic element
        N = 2*N;
        Pb = 0:1/N:1;
        Tb = [1:2:N-1; 3:2:N+1; 2:2:N];
    case 103   % cubic element
        N = 3*N;
        Pb = 0:1/N:1;
        Tb = [1:3:N-2; 4:3:N+1; 2:3:N-1; 3:3:N];
end
end

