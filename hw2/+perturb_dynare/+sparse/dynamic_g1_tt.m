function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = perturb_dynare.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 20
    T = [T; NaN(20 - size(T, 1), 1)];
end
T(10) = exp(y(8))*.67*(-1)/(y(5)*y(5));
T(11) = .97*(-1)/(y(9)*y(9));
T(12) = (-y(7))/(y(6)*y(6));
T(13) = .33*T(4)^(-0.6699999999999999);
T(14) = .67*y(6)^(-0.33);
T(15) = (-y(7))/(y(10)*y(10));
T(16) = getPowerDeriv(T(1),(-.67),1);
T(17) = exp(y(8))*.33*y(3)^(-0.6699999999999999);
T(18) = 1/y(10);
T(19) = .33*exp(y(12))*T(16)*T(18);
T(20) = 1/y(6);
end
