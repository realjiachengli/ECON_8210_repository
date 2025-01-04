function [T_order, T] = dynamic_g2_tt(y, x, params, steady_state, T_order, T)
if T_order >= 2
    return
end
[T_order, T] = perturb_dynare.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
T_order = 2;
if size(T, 1) < 36
    T = [T; NaN(36 - size(T, 1), 1)];
end
T(21) = .97*(y(9)+y(9))/(y(9)*y(9)*y(9)*y(9));
T(22) = (-((-y(7))*(y(10)+y(10))))/(y(10)*y(10)*y(10)*y(10));
T(23) = getPowerDeriv(T(1),(-.67),2);
T(24) = .33*exp(y(12))*(T(16)*T(22)+T(15)*T(15)*T(23));
T(25) = (-1)/(y(10)*y(10));
T(26) = .33*exp(y(12))*(T(16)*T(25)+T(15)*T(18)*T(23));
T(27) = .33*exp(y(12))*T(18)*T(18)*T(23);
T(28) = exp(y(8))*.67*(y(5)+y(5))/(y(5)*y(5)*y(5)*y(5));
T(29) = (-((-y(7))*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6));
T(30) = getPowerDeriv(T(4),(-0.6699999999999999),1);
T(31) = T(13)*T(29)+T(12)*.33*T(12)*T(30);
T(32) = (-1)/(y(6)*y(6));
T(33) = T(13)*T(32)+T(12)*.33*T(20)*T(30);
T(34) = T(20)*.33*T(20)*T(30);
T(35) = .67*getPowerDeriv(y(6),(-0.33),1);
T(36) = exp(y(8))*.33*getPowerDeriv(y(3),(-0.6699999999999999),1);
end
