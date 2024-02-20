function result = FE_basis_local_fun_1D(x,vertices,basis_type,basis_index,basis_der)

xl = vertices(1);     %left. 
xr = vertices(2);     %right.

h = xl - xr;

if basis_type == 101 % for linear element(P1).
    switch basis_index
        case 1
            switch basis_der
                case 0
                    result = (x - xr)/(xl - xr); % use lagrange form, insetead form from PDF(easy to extended).
                case 1
                    result = -1/h*ones(1,length(x)); % for simplicity,use vector form(in Gauss_quad_1d.m ).
                case 2
                    result = 0;
                otherwise
                    warning('All derivative >=2 is 0, confirm your input again. Use 2 for higher order.');
            end
        case 2
            switch basis_der
                case 0
                    result = (x - xl)/(xr - xl);
                case 1
                    result = 1/h*ones(1,length(x));
                case 2
                    result = 0;
                otherwise
                    warning('All derivative >=2 is 0, confirm your input again. Use 2 for higher order.');
            end
    end

elseif basis_type == 102 % for quadratic element(P2).
    xm = (xl+xr)/2;  %middle, quadratic element will use middle point.
    switch basis_index
        case 1  % for xl
            switch basis_der
                case 0
                    result = (x-xm).*(x-xr)/(xl-xm)/(xl-xr);
                case 1
                    result = (x-xr + x-xm)/(xl-xm)/(xl-xr);
                case 2
                    result = 2/(xl-xm)/(xl-xr)*ones(1,length(x));
                case 3
                    result = 0;
                otherwise
                    warning('All derivative >=3 is 0, confirm your input again. Use 3 for higher order.');
            end
        case 2  % for xr
            switch basis_der
                case 0
                    result = (x-xl).*(x-xm)/(xr-xl)/(xr-xm);
                case 1
                    result = (x-xl + x-xm)/(xr-xl)/(xr-xm);
                case 2
                    result = 2/(xr-xl)/(xr-xm)*ones(1,length(x));
                case 3
                    result = 0;
                otherwise
                    warning('All derivative >=3 is 0, confirm your input again. Use 3 for higher order.');
            end
        case 3  % for xm
            switch basis_der
                case 0
                    result = (x-xl).*(x-xr)/(xm-xl)/(xm-xr);
                case 1
                    result = (x-xl + x-xr)/(xm-xl)/(xm-xr);
                case 2
                    result = 2/(xm-xl)/(xm-xr)*ones(1,length(x));
                case 3
                    result = 0;
                otherwise
                    warning('All derivative >=3 is 0, confirm your input again. Use 3 for higher order.');
            end
    end

elseif basis_type == 103  % for cubic element(P3).
    xb = 2/3*xl + 1/3*xr;  % cubic element need other two point.
    xt = 2/3*xr + 1/3*xl;  % xb for x_bottom(1/3 bottom), xt for x_top(2/3top).
    switch basis_index
        case 1  % for xl
            switch basis_der
                case 0
                    result = (x-xb).*(x-xt).*(x-xr)/(xl-xb)/(xl-xt)/(xl-xr);
                case 1
                    result = ((x-xt).*(x-xr) + (x-xb).*(x-xr) + (x-xb).*(x-xt))/(xl-xb)/(xl-xt)/(xl-xr);
                case 2
                    result = ((x-xt)+(x-xr) + (x-xb)+(x-xr) + (x-xb)+(x-xt))/(xl-xb)/(xl-xt)/(xl-xr);
                case 3
                    result = 6/(xl-xb)/(xl-xt)/(xl-xr);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
        case 2  % for xr
            switch basis_der
                case 0
                    result = (x-xl).*(x-xb).*(x-xt)/(xr-xl)/(xr-xb)/(xr-xt);
                case 1
                    result = ((x-xl).*(x-xb) + (x-xb).*(x-xt) + (x-xl).*(x-xt))/(xr-xl)/(xr-xb)/(xr-xt);
                case 2
                    result = ((x-xl)+(x-xb) + (x-xb)+(x-xt) + (x-xl)+(x-xt))/(xr-xl)/(xr-xb)/(xr-xt);
                case 3
                    result = 6/(xr-xl)/(xr-xb)/(xr-xt);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
        case 3  % for xb
            switch basis_der
                case 0
                    result = (x-xl).*(x-xt).*(x-xr)/(xb-xl)/(xb-xt)/(xb-xr);
                case 1
                    result = ((x-xl).*(x-xt) + (x-xt).*(x-xr) + (x-xl).*(x-xr))/(xb-xl)/(xb-xt)/(xb-xr);
                case 2
                    result = ((x-xl)+(x-xt) + (x-xt)+(x-xr) + (x-xl)+(x-xr))/(xb-xl)/(xb-xt)/(xb-xr);
                case 3
                    result = 6/(xb-xl)/(xb-xt)/(xb-xr);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
        case 4  % for xt
            switch basis_der
                case 0
                    result = (x-xl).*(x-xb).*(x-xr)/(xt-xl)/(xt-xb)/(xt-xr);
                case 1
                    result = ((x-xl).*(x-xb) + (x-xb).*(x-xr) + (x-xl).*(x-xr))/(xt-xl)/(xt-xb)/(xt-xr);
                case 2
                    result = ((x-xl)+(x-xb) + (x-xb)+(x-xr) + (x-xl)+(x-xr))/(xt-xl)/(xt-xb)/(xt-xr);
                case 3
                    result = 6/(xt-xl)/(xt-xb)/(xt-xr);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
    end
else
    warining('I dont construct this kind of element.');
end 
end

