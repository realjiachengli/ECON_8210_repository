function [T_order, T] = dynamic_g3_tt(y, x, params, steady_state, T_order, T)
if T_order >= 3
    return
end
[T_order, T] = perturb_dynare.sparse.dynamic_g2_tt(y, x, params, steady_state, T_order, T);
T_order = 3;
if size(T, 1) < 40
    T = [T; NaN(40 - size(T, 1), 1)];
end
T(37) = (2*y(5)*y(5)*y(5)*y(5)-(y(5)+y(5))*(y(5)*y(5)*(y(5)+y(5))+y(5)*y(5)*(y(5)+y(5))))/(y(5)*y(5)*y(5)*y(5)*y(5)*y(5)*y(5)*y(5));
T(38) = getPowerDeriv(T(1),(-.67),3);
T(39) = T(29)*.33*T(12)*T(30);
T(40) = getPowerDeriv(T(4),(-0.6699999999999999),2);
end
