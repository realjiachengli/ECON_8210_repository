function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 36);

T = perturb_dynare.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(21) = .97*(y(7)+y(7))/(y(7)*y(7)*y(7)*y(7));
T(22) = (-((-y(5))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8));
T(23) = getPowerDeriv(T(1),(-.67),2);
T(24) = .33*exp(y(9))*(T(16)*T(22)+T(15)*T(15)*T(23));
T(25) = (-1)/(y(8)*y(8));
T(26) = .33*exp(y(9))*(T(16)*T(25)+T(15)*T(18)*T(23));
T(27) = .33*exp(y(9))*T(18)*T(18)*T(23);
T(28) = exp(y(6))*.67*(y(3)+y(3))/(y(3)*y(3)*y(3)*y(3));
T(29) = (-((-y(5))*(y(4)+y(4))))/(y(4)*y(4)*y(4)*y(4));
T(30) = getPowerDeriv(T(4),(-0.6699999999999999),1);
T(31) = T(13)*T(29)+T(12)*.33*T(12)*T(30);
T(32) = (-1)/(y(4)*y(4));
T(33) = T(13)*T(32)+T(12)*.33*T(20)*T(30);
T(34) = T(20)*.33*T(20)*T(30);
T(35) = .67*getPowerDeriv(y(4),(-0.33),1);
T(36) = exp(y(6))*.33*getPowerDeriv(y(1),(-0.6699999999999999),1);

end
