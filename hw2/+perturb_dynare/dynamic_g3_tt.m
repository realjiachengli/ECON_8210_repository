function T = dynamic_g3_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g3_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 40);

T = perturb_dynare.dynamic_g2_tt(T, y, x, params, steady_state, it_);

T(37) = (2*y(3)*y(3)*y(3)*y(3)-(y(3)+y(3))*(y(3)*y(3)*(y(3)+y(3))+y(3)*y(3)*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3)*y(3)*y(3)*y(3)*y(3));
T(38) = getPowerDeriv(T(1),(-.67),3);
T(39) = T(29)*.33*T(12)*T(30);
T(40) = getPowerDeriv(T(4),(-0.6699999999999999),2);

end
