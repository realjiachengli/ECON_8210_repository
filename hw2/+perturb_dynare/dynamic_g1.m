function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = perturb_dynare.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(4, 10);
g1(1,2)=(-.95);
g1(1,6)=1;
g1(1,10)=(-.007);
g1(2,3)=(-1)/(y(3)*y(3));
g1(2,7)=(-(T(3)*T(11)));
g1(2,8)=(-(T(8)*.33*exp(y(9))*T(15)*T(16)));
g1(2,5)=(-(T(8)*T(19)));
g1(2,9)=(-(T(2)*T(8)));
g1(3,3)=(-(T(5)*T(10)));
g1(3,4)=1-T(9)*T(12)*T(13);
g1(3,5)=(-(T(9)*T(13)*T(20)));
g1(3,6)=(-(T(5)*T(9)));
g1(4,3)=1;
g1(4,4)=(-(T(6)*T(14)));
g1(4,1)=(-(.9+T(7)*T(17)));
g1(4,5)=1;
g1(4,6)=(-(T(6)*T(7)));

end
