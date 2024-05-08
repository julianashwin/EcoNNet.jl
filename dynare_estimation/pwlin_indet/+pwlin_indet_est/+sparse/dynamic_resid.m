function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = pwlin_indet_est.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(5, 1);
    residual(1) = (y(7)) - (params(1)*y(8)+params(2)*y(6)+y(10));
    residual(2) = (y(7)-y(3)) - (x(3));
    residual(3) = (y(6)) - (params(3)*y(1)-params(4)*(y(7)*params(5)-y(8))+y(9));
    residual(4) = (y(10)) - (params(6)*y(5)+x(2));
    residual(5) = (y(9)) - (params(7)*y(4)+x(1));
end
