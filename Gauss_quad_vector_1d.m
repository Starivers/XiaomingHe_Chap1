function result = Gauss_quad_vector_1d(gauss_point, gauss_weight, vertices, function_name, test_basis_type, test_basis_index, test_basis_der)  
% The same as Gauss_quad_matrix_1d.m, except for without trial element.
xl = vertices(1);
xr = vertices(2);
% change gauss point from reference element to local element( from [-1,1] to [xl,xr] ).
gauss_point_local  = xl + (xr - xl)*(1+gauss_point)/2;
gauss_weight_local = gauss_weight*(xr - xl)/2;

test_val  = FE_basis_local_fun_1D(gauss_point_local, vertices, test_basis_type, test_basis_index, test_basis_der);
func_val  = feval(function_name, gauss_point_local);

result = gauss_weight_local*(test_val.*func_val)';
end
