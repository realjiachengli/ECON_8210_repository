function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 8);

T(1) = y(3)/y(2);
T(2) = .33*exp(y(4))*T(1)^(-.67);
T(3) = T(2)+.9;
T(4) = 1/y(1)*.97;
T(5) = T(1)^.33;
T(6) = exp(y(4))*1/y(1)*.67;
T(7) = exp(y(4))*y(3)^.33;
T(8) = y(2)^.67;

end
