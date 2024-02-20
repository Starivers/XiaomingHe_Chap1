function err = FE_solver_1d_possion(num_of_element,gauss_type,trial_basis_type,test_basis_type)
% FE_solver_1d_possion(4,3,101,101)
% FE_solver_1d_possion(4,3,102,102)
% FE_solver_1d_possion(4,4,103,103) % here 103 is cubic element, instead of Hermite element in Pdf.
% This is a code for lecture of Pro Xiaoming He.
% The detailed equation refer to the Pdf "Chapter 1: Finite elements for 1D second order elliptic equation"

% meaning of input---------------------
% num_of_element: number of element.
% gauss_type: how many gauss points to use within a integral.
% basis_type: which basis to use. 101 for linear, 102 for quadratic, 103 for cubic.
% -------------------------------------
[P,T] = generate_PbTb(num_of_element,101);

[P_trial,T_trial] = generate_PbTb(num_of_element, trial_basis_type);
[P_test,  T_test] = generate_PbTb(num_of_element,  test_basis_type);

matrix_size = [length(P_test),length(P_trial)];
vector_size = length(P_trial);

[gauss_weight,gauss_point] = gaussValues_1d(gauss_type); % both gauss_point and weight are row vector.

trial_basis_der_A = 1;
test_basis_der_A  = 1;
A = assemble_matrix_1d(matrix_size, gauss_point, gauss_weight, 'function_c', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A, test_basis_type, test_basis_der_A);    

test_basis_der_b = 0;
b = assemble_vector_1d(vector_size, gauss_point, gauss_weight, 'function_f', P, T, T_test, test_basis_type, test_basis_der_b);

[A,b] = boundary_process(A,b);

sol = A\b;

plot(P_trial,sol,'o','lineWidth',1)
hold on 
plot(P_trial,P_trial.*cos(P_trial),'-*','lineWidth',1)
%err = norm(sol'-P_trial.*cos(P_trial))/sqrt(num_of_element);  % L2 error
err = max(sol'-P_trial.*cos(P_trial));                         % L_inf error
legend("numerical sol","real sol")
end


