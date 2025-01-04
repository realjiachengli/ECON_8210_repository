function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 9
    T = [T; NaN(9 - size(T, 1), 1)];
end
T(1) = y(7)/y(10);
T(2) = .33*exp(y(12))*T(1)^(-.67);
T(3) = T(2)+.9;
T(4) = y(7)/y(6);
T(5) = T(4)^.33;
T(6) = exp(y(8))*y(3)^.33;
T(7) = y(6)^.67;
T(8) = 1/y(9)*.97;
T(9) = exp(y(8))*1/y(5)*.67;
end
