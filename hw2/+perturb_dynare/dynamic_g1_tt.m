function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 20);

T = perturb_dynare.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(10) = exp(y(6))*.67*(-1)/(y(3)*y(3));
T(11) = .97*(-1)/(y(7)*y(7));
T(12) = (-y(5))/(y(4)*y(4));
T(13) = .33*T(4)^(-0.6699999999999999);
T(14) = .67*y(4)^(-0.33);
T(15) = (-y(5))/(y(8)*y(8));
T(16) = getPowerDeriv(T(1),(-.67),1);
T(17) = exp(y(6))*.33*y(1)^(-0.6699999999999999);
T(18) = 1/y(8);
T(19) = .33*exp(y(9))*T(16)*T(18);
T(20) = 1/y(4);

end
