function [y, T, residual, g1] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(3, 1);
  T(1)=exp(y(4));
  T(2)=.33*T(1)*(y(3)/y(2))^(-.67)+.9;
  residual(1)=(1/y(1))-(T(2)*1/y(1)*.97);
  T(3)=(y(3)/y(2))^.33;
  residual(2)=(y(2))-(T(3)*T(1)*1/y(1)*.67);
  T(4)=T(1)*y(3)^.33;
  T(5)=y(2)^.67;
  residual(3)=(y(3))-(T(4)*T(5)+y(3)*.9-y(1));
  T(6)=getPowerDeriv(y(3)/y(2),(-.67),1);
  T(7)=.33*(y(3)/y(2))^(-0.6699999999999999);
if nargout > 3
    g1_v = NaN(9, 1);
g1_v(1)=(-(1/y(1)*.97*.33*T(1)*T(6)*1/y(2)));
g1_v(2)=(-(T(1)*1/y(1)*.67*T(7)*1/y(2)));
g1_v(3)=1-(.9+T(5)*T(1)*.33*y(3)^(-0.6699999999999999));
g1_v(4)=(-(1/y(1)*.97*.33*T(1)*(-y(3))/(y(2)*y(2))*T(6)));
g1_v(5)=1-T(1)*1/y(1)*.67*(-y(3))/(y(2)*y(2))*T(7);
g1_v(6)=(-(T(4)*.67*y(2)^(-0.33)));
g1_v(7)=(-1)/(y(1)*y(1))-T(2)*.97*(-1)/(y(1)*y(1));
g1_v(8)=(-(T(3)*T(1)*.67*(-1)/(y(1)*y(1))));
g1_v(9)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 3, 3);
end
end
