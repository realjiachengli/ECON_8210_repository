function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(8)=.95*y(4)+.007*x(1);
end
