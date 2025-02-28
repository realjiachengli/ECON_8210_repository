function [g3_v, T_order, T] = dynamic_g3(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(40, 1);
end
[T_order, T] = perturb_dynare.sparse.dynamic_g3_tt(y, x, params, steady_state, T_order, T);
g3_v = NaN(51, 1);
g3_v(1)=T(37);
g3_v(2)=(-(T(3)*.97*(2*y(9)*y(9)*y(9)*y(9)-(y(9)+y(9))*(y(9)*y(9)*(y(9)+y(9))+y(9)*y(9)*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9)*y(9)*y(9)*y(9)*y(9))));
g3_v(3)=(-(.33*exp(y(12))*T(15)*T(16)*T(21)));
g3_v(4)=(-(T(19)*T(21)));
g3_v(5)=(-(T(2)*T(21)));
g3_v(6)=(-(T(11)*T(24)));
g3_v(7)=(-(T(11)*T(26)));
g3_v(8)=(-(T(11)*.33*exp(y(12))*T(15)*T(16)));
g3_v(9)=(-(T(11)*T(27)));
g3_v(10)=(-(T(11)*T(19)));
g3_v(11)=(-(T(2)*T(11)));
g3_v(12)=(-(T(8)*.33*exp(y(12))*(T(22)*T(15)*T(23)+T(16)*(y(10)*y(10)*y(10)*y(10)*(-(2*(-y(7))))-(-((-y(7))*(y(10)+y(10))))*(y(10)*y(10)*(y(10)+y(10))+y(10)*y(10)*(y(10)+y(10))))/(y(10)*y(10)*y(10)*y(10)*y(10)*y(10)*y(10)*y(10))+T(22)*T(15)*T(23)+T(15)*(T(22)*T(23)+T(15)*T(15)*T(38)))));
g3_v(13)=(-(T(8)*.33*exp(y(12))*(T(22)*T(18)*T(23)+T(16)*(y(10)+y(10))/(y(10)*y(10)*y(10)*y(10))+T(15)*T(23)*T(25)+T(15)*(T(23)*T(25)+T(15)*T(18)*T(38)))));
g3_v(14)=(-(T(8)*T(24)));
g3_v(15)=(-(T(8)*.33*exp(y(12))*(T(25)*T(18)*T(23)+T(25)*T(18)*T(23)+T(15)*T(18)*T(18)*T(38))));
g3_v(16)=(-(T(8)*T(26)));
g3_v(17)=(-(T(8)*.33*exp(y(12))*T(15)*T(16)));
g3_v(18)=(-(T(8)*.33*exp(y(12))*T(18)*T(18)*T(18)*T(38)));
g3_v(19)=(-(T(8)*T(27)));
g3_v(20)=(-(T(8)*T(19)));
g3_v(21)=(-(T(2)*T(8)));
g3_v(22)=(-(T(5)*exp(y(8))*.67*T(37)));
g3_v(23)=(-(T(12)*T(13)*T(28)));
g3_v(24)=(-(T(13)*T(20)*T(28)));
g3_v(25)=(-(T(5)*T(28)));
g3_v(26)=(-(T(10)*T(31)));
g3_v(27)=(-(T(10)*T(33)));
g3_v(28)=(-(T(10)*T(12)*T(13)));
g3_v(29)=(-(T(10)*T(34)));
g3_v(30)=(-(T(10)*T(13)*T(20)));
g3_v(31)=(-(T(5)*T(10)));
g3_v(32)=(-(T(9)*(T(39)+T(13)*(y(6)*y(6)*y(6)*y(6)*(-(2*(-y(7))))-(-((-y(7))*(y(6)+y(6))))*(y(6)*y(6)*(y(6)+y(6))+y(6)*y(6)*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6)*y(6)*y(6)*y(6)*y(6))+T(39)+T(12)*.33*(T(29)*T(30)+T(12)*T(12)*T(40)))));
g3_v(33)=(-(T(9)*(T(29)*.33*T(20)*T(30)+T(13)*(y(6)+y(6))/(y(6)*y(6)*y(6)*y(6))+.33*T(12)*T(30)*T(32)+T(12)*.33*(T(30)*T(32)+T(12)*T(20)*T(40)))));
g3_v(34)=(-(T(9)*T(31)));
g3_v(35)=(-(T(9)*(T(32)*.33*T(20)*T(30)+T(32)*.33*T(20)*T(30)+T(12)*.33*T(20)*T(20)*T(40))));
g3_v(36)=(-(T(9)*T(33)));
g3_v(37)=(-(T(9)*T(12)*T(13)));
g3_v(38)=(-(T(9)*T(20)*.33*T(20)*T(20)*T(40)));
g3_v(39)=(-(T(9)*T(34)));
g3_v(40)=(-(T(9)*T(13)*T(20)));
g3_v(41)=(-(T(5)*T(9)));
g3_v(42)=(-(T(6)*.67*getPowerDeriv(y(6),(-0.33),2)));
g3_v(43)=(-(T(17)*T(35)));
g3_v(44)=(-(T(6)*T(35)));
g3_v(45)=(-(T(14)*T(36)));
g3_v(46)=(-(T(14)*T(17)));
g3_v(47)=(-(T(6)*T(14)));
g3_v(48)=(-(T(7)*exp(y(8))*.33*getPowerDeriv(y(3),(-0.6699999999999999),2)));
g3_v(49)=(-(T(7)*T(36)));
g3_v(50)=(-(T(7)*T(17)));
g3_v(51)=(-(T(6)*T(7)));
end
