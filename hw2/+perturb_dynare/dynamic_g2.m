function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = perturb_dynare.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(42,1);
g2_j = zeros(42,1);
g2_v = zeros(42,1);

g2_i(1)=2;
g2_i(2)=2;
g2_i(3)=2;
g2_i(4)=2;
g2_i(5)=2;
g2_i(6)=2;
g2_i(7)=2;
g2_i(8)=2;
g2_i(9)=2;
g2_i(10)=2;
g2_i(11)=2;
g2_i(12)=2;
g2_i(13)=2;
g2_i(14)=2;
g2_i(15)=2;
g2_i(16)=2;
g2_i(17)=2;
g2_i(18)=3;
g2_i(19)=3;
g2_i(20)=3;
g2_i(21)=3;
g2_i(22)=3;
g2_i(23)=3;
g2_i(24)=3;
g2_i(25)=3;
g2_i(26)=3;
g2_i(27)=3;
g2_i(28)=3;
g2_i(29)=3;
g2_i(30)=3;
g2_i(31)=3;
g2_i(32)=3;
g2_i(33)=3;
g2_i(34)=4;
g2_i(35)=4;
g2_i(36)=4;
g2_i(37)=4;
g2_i(38)=4;
g2_i(39)=4;
g2_i(40)=4;
g2_i(41)=4;
g2_i(42)=4;
g2_j(1)=23;
g2_j(2)=67;
g2_j(3)=68;
g2_j(4)=77;
g2_j(5)=65;
g2_j(6)=47;
g2_j(7)=69;
g2_j(8)=87;
g2_j(9)=78;
g2_j(10)=75;
g2_j(11)=48;
g2_j(12)=79;
g2_j(13)=88;
g2_j(14)=45;
g2_j(15)=49;
g2_j(16)=85;
g2_j(17)=89;
g2_j(18)=23;
g2_j(19)=24;
g2_j(20)=33;
g2_j(21)=25;
g2_j(22)=43;
g2_j(23)=26;
g2_j(24)=53;
g2_j(25)=34;
g2_j(26)=35;
g2_j(27)=44;
g2_j(28)=36;
g2_j(29)=54;
g2_j(30)=45;
g2_j(31)=46;
g2_j(32)=55;
g2_j(33)=56;
g2_j(34)=34;
g2_j(35)=31;
g2_j(36)=4;
g2_j(37)=36;
g2_j(38)=54;
g2_j(39)=1;
g2_j(40)=6;
g2_j(41)=51;
g2_j(42)=56;
g2_v(1)=(y(3)+y(3))/(y(3)*y(3)*y(3)*y(3));
g2_v(2)=(-(T(3)*T(21)));
g2_v(3)=(-(T(11)*.33*exp(y(9))*T(15)*T(16)));
g2_v(4)=g2_v(3);
g2_v(5)=(-(T(11)*T(19)));
g2_v(6)=g2_v(5);
g2_v(7)=(-(T(2)*T(11)));
g2_v(8)=g2_v(7);
g2_v(9)=(-(T(8)*T(24)));
g2_v(10)=(-(T(8)*T(26)));
g2_v(11)=g2_v(10);
g2_v(12)=(-(T(8)*.33*exp(y(9))*T(15)*T(16)));
g2_v(13)=g2_v(12);
g2_v(14)=(-(T(8)*T(27)));
g2_v(15)=(-(T(8)*T(19)));
g2_v(16)=g2_v(15);
g2_v(17)=(-(T(2)*T(8)));
g2_v(18)=(-(T(5)*T(28)));
g2_v(19)=(-(T(10)*T(12)*T(13)));
g2_v(20)=g2_v(19);
g2_v(21)=(-(T(10)*T(13)*T(20)));
g2_v(22)=g2_v(21);
g2_v(23)=(-(T(5)*T(10)));
g2_v(24)=g2_v(23);
g2_v(25)=(-(T(9)*T(31)));
g2_v(26)=(-(T(9)*T(33)));
g2_v(27)=g2_v(26);
g2_v(28)=(-(T(9)*T(12)*T(13)));
g2_v(29)=g2_v(28);
g2_v(30)=(-(T(9)*T(34)));
g2_v(31)=(-(T(9)*T(13)*T(20)));
g2_v(32)=g2_v(31);
g2_v(33)=(-(T(5)*T(9)));
g2_v(34)=(-(T(6)*T(35)));
g2_v(35)=(-(T(14)*T(17)));
g2_v(36)=g2_v(35);
g2_v(37)=(-(T(6)*T(14)));
g2_v(38)=g2_v(37);
g2_v(39)=(-(T(7)*T(36)));
g2_v(40)=(-(T(7)*T(17)));
g2_v(41)=g2_v(40);
g2_v(42)=(-(T(6)*T(7)));
g2 = sparse(g2_i,g2_j,g2_v,4,100);
end
