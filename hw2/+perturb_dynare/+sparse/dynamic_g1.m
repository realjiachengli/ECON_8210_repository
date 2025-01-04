function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(20, 1);
end
[T_order, T] = perturb_dynare.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(17, 1);
g1_v(1)=(-(.9+T(7)*T(17)));
g1_v(2)=(-.95);
g1_v(3)=(-1)/(y(5)*y(5));
g1_v(4)=(-(T(5)*T(10)));
g1_v(5)=1;
g1_v(6)=1-T(9)*T(12)*T(13);
g1_v(7)=(-(T(6)*T(14)));
g1_v(8)=(-(T(8)*T(19)));
g1_v(9)=(-(T(9)*T(13)*T(20)));
g1_v(10)=1;
g1_v(11)=1;
g1_v(12)=(-(T(5)*T(9)));
g1_v(13)=(-(T(6)*T(7)));
g1_v(14)=(-(T(3)*T(11)));
g1_v(15)=(-(T(8)*.33*exp(y(12))*T(15)*T(16)));
g1_v(16)=(-(T(2)*T(8)));
g1_v(17)=(-.007);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 4, 13);
end
