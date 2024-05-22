function [y, T, residual, g1] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  residual(1)=(y(2))-(y(2)*params(1)+params(2)*y(1)+y(4));
  residual(2)=(y(1))-(y(1)*params(3)-params(4)*(y(2)*params(5)-y(2))+y(3));
if nargout > 3
    g1_v = NaN(4, 1);
g1_v(1)=1-params(1);
g1_v(2)=params(4)*(params(5)-1);
g1_v(3)=(-params(2));
g1_v(4)=1-params(3);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 2);
end
end
