function g3 = dynamic_g3(T, y, x, params, steady_state, it_, T_flag)
% function g3 = dynamic_g3(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g3
%

if T_flag
    T = perturb_dynare.dynamic_g3_tt(T, y, x, params, steady_state, it_);
end
g3_i = zeros(51,1);
g3_j = zeros(51,1);
g3_v = zeros(51,1);

g3_i(1)=2;
g3_i(2)=2;
g3_i(3)=2;
g3_i(4)=2;
g3_i(5)=2;
g3_i(6)=2;
g3_i(7)=2;
g3_i(8)=2;
g3_i(9)=2;
g3_i(10)=2;
g3_i(11)=2;
g3_i(12)=2;
g3_i(13)=2;
g3_i(14)=2;
g3_i(15)=2;
g3_i(16)=2;
g3_i(17)=2;
g3_i(18)=2;
g3_i(19)=2;
g3_i(20)=2;
g3_i(21)=2;
g3_i(22)=3;
g3_i(23)=3;
g3_i(24)=3;
g3_i(25)=3;
g3_i(26)=3;
g3_i(27)=3;
g3_i(28)=3;
g3_i(29)=3;
g3_i(30)=3;
g3_i(31)=3;
g3_i(32)=3;
g3_i(33)=3;
g3_i(34)=3;
g3_i(35)=3;
g3_i(36)=3;
g3_i(37)=3;
g3_i(38)=3;
g3_i(39)=3;
g3_i(40)=3;
g3_i(41)=3;
g3_i(42)=4;
g3_i(43)=4;
g3_i(44)=4;
g3_i(45)=4;
g3_i(46)=4;
g3_i(47)=4;
g3_i(48)=4;
g3_i(49)=4;
g3_i(50)=4;
g3_i(51)=4;
g3_j(1)=223;
g3_j(2)=667;
g3_j(3)=668;
g3_j(4)=665;
g3_j(5)=669;
g3_j(6)=678;
g3_j(7)=675;
g3_j(8)=679;
g3_j(9)=645;
g3_j(10)=649;
g3_j(11)=689;
g3_j(12)=778;
g3_j(13)=775;
g3_j(14)=779;
g3_j(15)=745;
g3_j(16)=749;
g3_j(17)=789;
g3_j(18)=445;
g3_j(19)=449;
g3_j(20)=489;
g3_j(21)=889;
g3_j(22)=223;
g3_j(23)=224;
g3_j(24)=225;
g3_j(25)=226;
g3_j(26)=234;
g3_j(27)=235;
g3_j(28)=236;
g3_j(29)=245;
g3_j(30)=246;
g3_j(31)=256;
g3_j(32)=334;
g3_j(33)=335;
g3_j(34)=336;
g3_j(35)=345;
g3_j(36)=346;
g3_j(37)=356;
g3_j(38)=445;
g3_j(39)=446;
g3_j(40)=456;
g3_j(41)=556;
g3_j(42)=334;
g3_j(43)=331;
g3_j(44)=336;
g3_j(45)=301;
g3_j(46)=306;
g3_j(47)=356;
g3_j(48)=1;
g3_j(49)=6;
g3_j(50)=56;
g3_j(51)=556;
g3_v(1)=T(37);
g3_v(2)=(-(T(3)*.97*(2*y(7)*y(7)*y(7)*y(7)-(y(7)+y(7))*(y(7)*y(7)*(y(7)+y(7))+y(7)*y(7)*(y(7)+y(7))))/(y(7)*y(7)*y(7)*y(7)*y(7)*y(7)*y(7)*y(7))));
g3_v(3)=(-(.33*exp(y(9))*T(15)*T(16)*T(21)));
g3_v(4)=(-(T(19)*T(21)));
g3_v(5)=(-(T(2)*T(21)));
g3_v(6)=(-(T(11)*T(24)));
g3_v(7)=(-(T(11)*T(26)));
g3_v(8)=(-(T(11)*.33*exp(y(9))*T(15)*T(16)));
g3_v(9)=(-(T(11)*T(27)));
g3_v(10)=(-(T(11)*T(19)));
g3_v(11)=(-(T(2)*T(11)));
g3_v(12)=(-(T(8)*.33*exp(y(9))*(T(22)*T(15)*T(23)+T(16)*(y(8)*y(8)*y(8)*y(8)*(-(2*(-y(5))))-(-((-y(5))*(y(8)+y(8))))*(y(8)*y(8)*(y(8)+y(8))+y(8)*y(8)*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8))+T(22)*T(15)*T(23)+T(15)*(T(22)*T(23)+T(15)*T(15)*T(38)))));
g3_v(13)=(-(T(8)*.33*exp(y(9))*(T(22)*T(18)*T(23)+T(16)*(y(8)+y(8))/(y(8)*y(8)*y(8)*y(8))+T(15)*T(23)*T(25)+T(15)*(T(23)*T(25)+T(15)*T(18)*T(38)))));
g3_v(14)=(-(T(8)*T(24)));
g3_v(15)=(-(T(8)*.33*exp(y(9))*(T(25)*T(18)*T(23)+T(25)*T(18)*T(23)+T(15)*T(18)*T(18)*T(38))));
g3_v(16)=(-(T(8)*T(26)));
g3_v(17)=(-(T(8)*.33*exp(y(9))*T(15)*T(16)));
g3_v(18)=(-(T(8)*.33*exp(y(9))*T(18)*T(18)*T(18)*T(38)));
g3_v(19)=(-(T(8)*T(27)));
g3_v(20)=(-(T(8)*T(19)));
g3_v(21)=(-(T(2)*T(8)));
g3_v(22)=(-(T(5)*exp(y(6))*.67*T(37)));
g3_v(23)=(-(T(12)*T(13)*T(28)));
g3_v(24)=(-(T(13)*T(20)*T(28)));
g3_v(25)=(-(T(5)*T(28)));
g3_v(26)=(-(T(10)*T(31)));
g3_v(27)=(-(T(10)*T(33)));
g3_v(28)=(-(T(10)*T(12)*T(13)));
g3_v(29)=(-(T(10)*T(34)));
g3_v(30)=(-(T(10)*T(13)*T(20)));
g3_v(31)=(-(T(5)*T(10)));
g3_v(32)=(-(T(9)*(T(39)+T(13)*(y(4)*y(4)*y(4)*y(4)*(-(2*(-y(5))))-(-((-y(5))*(y(4)+y(4))))*(y(4)*y(4)*(y(4)+y(4))+y(4)*y(4)*(y(4)+y(4))))/(y(4)*y(4)*y(4)*y(4)*y(4)*y(4)*y(4)*y(4))+T(39)+T(12)*.33*(T(29)*T(30)+T(12)*T(12)*T(40)))));
g3_v(33)=(-(T(9)*(T(29)*.33*T(20)*T(30)+T(13)*(y(4)+y(4))/(y(4)*y(4)*y(4)*y(4))+.33*T(12)*T(30)*T(32)+T(12)*.33*(T(30)*T(32)+T(12)*T(20)*T(40)))));
g3_v(34)=(-(T(9)*T(31)));
g3_v(35)=(-(T(9)*(T(32)*.33*T(20)*T(30)+T(32)*.33*T(20)*T(30)+T(12)*.33*T(20)*T(20)*T(40))));
g3_v(36)=(-(T(9)*T(33)));
g3_v(37)=(-(T(9)*T(12)*T(13)));
g3_v(38)=(-(T(9)*T(20)*.33*T(20)*T(20)*T(40)));
g3_v(39)=(-(T(9)*T(34)));
g3_v(40)=(-(T(9)*T(13)*T(20)));
g3_v(41)=(-(T(5)*T(9)));
g3_v(42)=(-(T(6)*.67*getPowerDeriv(y(4),(-0.33),2)));
g3_v(43)=(-(T(17)*T(35)));
g3_v(44)=(-(T(6)*T(35)));
g3_v(45)=(-(T(14)*T(36)));
g3_v(46)=(-(T(14)*T(17)));
g3_v(47)=(-(T(6)*T(14)));
g3_v(48)=(-(T(7)*exp(y(6))*.33*getPowerDeriv(y(1),(-0.6699999999999999),2)));
g3_v(49)=(-(T(7)*T(36)));
g3_v(50)=(-(T(7)*T(17)));
g3_v(51)=(-(T(6)*T(7)));
g3 = sparse(g3_i,g3_j,g3_v,4,1000);
end
