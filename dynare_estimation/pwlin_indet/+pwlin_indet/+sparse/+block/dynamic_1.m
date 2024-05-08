function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(7)=y(3)+x(3);
  y(10)=params(6)*y(5)+x(2);
  y(9)=params(7)*y(4)+x(1);
end
