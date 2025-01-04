function [g1, T_order, T] = static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 8
    T_order = -1;
    T = NaN(10, 1);
end
[T_order, T] = perturb_dynare.sparse.static_g1_tt(y, x, params, T_order, T);
g1_v = NaN(13, 1);
g1_v(1)=(-1)/(y(1)*y(1))-T(3)*.97*(-1)/(y(1)*y(1));
g1_v(2)=(-(T(5)*exp(y(4))*.67*(-1)/(y(1)*y(1))));
g1_v(3)=1;
g1_v(4)=(-(T(4)*.33*exp(y(4))*(-y(3))/(y(2)*y(2))*T(9)));
g1_v(5)=1-T(6)*(-y(3))/(y(2)*y(2))*T(10);
g1_v(6)=(-(T(7)*.67*y(2)^(-0.33)));
g1_v(7)=(-(T(4)*.33*exp(y(4))*T(9)*1/y(2)));
g1_v(8)=(-(T(6)*T(10)*1/y(2)));
g1_v(9)=1-(.9+T(8)*exp(y(4))*.33*y(3)^(-0.6699999999999999));
g1_v(10)=0.05000000000000004;
g1_v(11)=(-(T(2)*T(4)));
g1_v(12)=(-(T(5)*T(6)));
g1_v(13)=(-(T(7)*T(8)));
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 4, 4);
end
