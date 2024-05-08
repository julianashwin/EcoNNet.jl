function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = pwlin_det_est_learn.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(15, 1);
g1_v(1)=(-params(3));
g1_v(2)=(-params(9));
g1_v(3)=(-params(8));
g1_v(4)=(-params(2));
g1_v(5)=1;
g1_v(6)=1;
g1_v(7)=params(4)*params(5);
g1_v(8)=(-1);
g1_v(9)=1;
g1_v(10)=(-1);
g1_v(11)=1;
g1_v(12)=(-params(1));
g1_v(13)=(-params(4));
g1_v(14)=(-1);
g1_v(15)=(-1);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 4, 14);
end
