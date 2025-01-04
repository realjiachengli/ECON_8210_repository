function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(3, 1);
  T(1)=exp(y(8));
  T(2)=T(1)*y(3)^.33;
  T(3)=y(6)^.67;
  residual(1)=(y(7))-(T(2)*T(3)+.9*y(3)-y(5));
  T(4)=.33*exp(y(12));
  T(5)=T(4)*(y(7)/y(10))^(-.67)+.9;
  residual(2)=(1/y(5))-(T(5)*1/y(9)*.97);
  T(6)=(y(7)/y(6))^.33;
  residual(3)=(y(6))-(T(6)*T(1)*1/y(5)*.67);
  T(7)=.33*(y(7)/y(6))^(-0.6699999999999999);
  T(8)=getPowerDeriv(y(7)/y(10),(-.67),1);
if nargout > 3
    g1_v = NaN(11, 1);
g1_v(1)=(-(.9+T(3)*T(1)*.33*y(3)^(-0.6699999999999999)));
g1_v(2)=1;
g1_v(3)=(-(1/y(9)*.97*T(4)*T(8)*1/y(10)));
g1_v(4)=(-(T(1)*1/y(5)*.67*T(7)*1/y(6)));
g1_v(5)=1;
g1_v(6)=(-1)/(y(5)*y(5));
g1_v(7)=(-(T(6)*T(1)*.67*(-1)/(y(5)*y(5))));
g1_v(8)=(-(T(2)*.67*y(6)^(-0.33)));
g1_v(9)=1-T(1)*1/y(5)*.67*(-y(7))/(y(6)*y(6))*T(7);
g1_v(10)=(-(T(5)*.97*(-1)/(y(9)*y(9))));
g1_v(11)=(-(1/y(9)*.97*T(4)*(-y(7))/(y(10)*y(10))*T(8)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 3, 9);
end
end
