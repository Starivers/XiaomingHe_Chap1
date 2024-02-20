function b = assemble_vector_1d(vector_size, gauss_point, gauss_weight, function_name, P, T, T_test, test_basis_type, test_basis_der)

b = sparse(vector_size,1);
num_local_test = size(T_test,1);

for n = 1:size(T,2)
    vertices = P(T(:,n));
    for beta = 1:num_local_test
        val = Gauss_quad_vector_1d(gauss_point, gauss_weight, vertices, function_name, test_basis_type, beta, test_basis_der);             
        b(T_test(beta,n)) = b(T_test(beta,n)) + val;
    end
end

end
