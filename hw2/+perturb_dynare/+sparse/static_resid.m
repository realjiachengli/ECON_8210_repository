function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = perturb_dynare.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(4, 1);
    residual(1) = (y(4)) - (y(4)*.95+.007*x(1));
    residual(2) = (1/y(1)) - (T(3)*T(4));
    residual(3) = (y(2)) - (T(5)*T(6));
    residual(4) = (y(3)) - (T(7)*T(8)+y(3)*.9-y(1));
end
