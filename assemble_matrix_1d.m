function A = assemble_matrix_1d(matrix_size, gauss_point, gauss_weight, function_name, P, T, T_trial, T_test, trial_basis_type, trial_basis_der, test_basis_type, test_basis_der)
% column of A corresponds to trial space, row of A corresponds to test space.

A = sparse(matrix_size(1),matrix_size(2));

num_local_trial = size(T_trial,1); % how many basis within an element.
num_local_test  = size(T_test ,1);

for n = 1:size(T,2)
    vertices = P(T(:,n));    % row vector. for higher dimension, use vertices = P(:,T(:,n)) instead(every column composites a coordinate).
    for alpha = 1:num_local_trial  % trial, corresponds to column.
        for beta = 1:num_local_test % test, correspinds to row.
            %Gauss_quad_1d(gauss_point,gauss_weight,vertices,'function_name',trial_basis_type,trial_basis_index,trial_basis_der, test_basis_type,test_basis_index,test_basis_der)
            val = Gauss_quad_matrix_1d(gauss_point, gauss_weight, vertices, function_name, trial_basis_type,alpha, trial_basis_der, test_basis_type,beta, test_basis_der);
            A(T_test(beta,n),T_trial(alpha,n)) = A(T_test(beta,n),T_trial(alpha,n)) + val;
        end
    end
end

end