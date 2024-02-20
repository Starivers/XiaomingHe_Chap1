function result = Gauss_quad_matrix_1d(gauss_point,gauss_weight,vertices,function_name,trial_basis_type,trial_basis_index,trial_basis_der, test_basis_type,test_basis_index,test_basis_der)
% This part is different from Pro Xiaoming He.
% He generates gauss point in local element. I use global reference gauss point, and use change of variable.
% In other word, my gauss_point and weight is global.
% This part you can refer to "https://numfactory.upc.edu/web/Calculo2/P2_Integracio/html/Gauss1D.html#H_954FFA45" for more detail.

xl = vertices(1);
xr = vertices(2);
% change gauss point from reference element to local element( from [-1,1] to [xl,xr] ).
gauss_point_local  = xl + (xr - xl)*(1+gauss_point)/2;
gauss_weight_local = gauss_weight*(xr - xl)/2;

trial_val = FE_basis_local_fun_1D(gauss_point_local,vertices,trial_basis_type,trial_basis_index,trial_basis_der);
test_val  = FE_basis_local_fun_1D(gauss_point_local,vertices, test_basis_type, test_basis_index, test_basis_der);
 
func_val  = feval(function_name,gauss_point_local);
% Code in bilibili: @function_name(gauss_point_local),which is wrong. How to use function handle @ to rewrite this part?

% input gauss_point and return of FE_basis_local_fun_1D are row vector.
result = gauss_weight_local*(trial_val.*test_val.*func_val)';
end

