function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(1) = y(3)/y(2);
T(2) = .33*exp(y(4))*T(1)^(-.67);
T(3) = T(2)+.9;
T(4) = 1/y(1)*.97;
T(5) = T(1)^.33;
T(6) = exp(y(4))*1/y(1)*.67;
T(7) = exp(y(4))*y(3)^.33;
T(8) = y(2)^.67;
end
