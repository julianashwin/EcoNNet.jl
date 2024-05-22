function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(8)=params(8)*y(4)+x(2);
  y(7)=params(9)*y(3)+x(1);
end
