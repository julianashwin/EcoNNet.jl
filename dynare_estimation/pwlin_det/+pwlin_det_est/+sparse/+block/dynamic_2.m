function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  residual(1)=(y(5))-(params(3)*y(1)-params(4)*(y(6)*params(5)-y(10))+y(7));
  residual(2)=(y(6))-(params(1)*y(10)+params(2)*y(5)+y(8));
if nargout > 3
    g1_v = NaN(7, 1);
g1_v(1)=(-params(3));
g1_v(2)=1;
g1_v(3)=(-params(2));
g1_v(4)=params(4)*params(5);
g1_v(5)=1;
g1_v(6)=(-params(4));
g1_v(7)=(-params(1));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 6);
end
end
