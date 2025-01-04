function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(9, 1);
end
[T_order, T] = perturb_dynare.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(4, 1);
    residual(1) = (y(8)) - (.95*y(4)+.007*x(1));
    residual(2) = (1/y(5)) - (T(3)*T(8));
    residual(3) = (y(6)) - (T(5)*T(9));
    residual(4) = (y(7)) - (T(6)*T(7)+.9*y(3)-y(5));
end
