function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = perturb_dynare.static_g1_tt(T, y, x, params);
end
g1 = zeros(4, 4);
g1(1,4)=0.05000000000000004;
g1(2,1)=(-1)/(y(1)*y(1))-T(3)*.97*(-1)/(y(1)*y(1));
g1(2,2)=(-(T(4)*.33*exp(y(4))*(-y(3))/(y(2)*y(2))*T(9)));
g1(2,3)=(-(T(4)*.33*exp(y(4))*T(9)*1/y(2)));
g1(2,4)=(-(T(2)*T(4)));
g1(3,1)=(-(T(5)*exp(y(4))*.67*(-1)/(y(1)*y(1))));
g1(3,2)=1-T(6)*(-y(3))/(y(2)*y(2))*T(10);
g1(3,3)=(-(T(6)*T(10)*1/y(2)));
g1(3,4)=(-(T(5)*T(6)));
g1(4,1)=1;
g1(4,2)=(-(T(7)*.67*y(2)^(-0.33)));
g1(4,3)=1-(.9+T(8)*exp(y(4))*.33*y(3)^(-0.6699999999999999));
g1(4,4)=(-(T(7)*T(8)));

end
