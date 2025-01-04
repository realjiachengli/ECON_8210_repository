function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = perturb_dynare.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 10
    T = [T; NaN(10 - size(T, 1), 1)];
end
T(9) = getPowerDeriv(T(1),(-.67),1);
T(10) = .33*T(1)^(-0.6699999999999999);
end
