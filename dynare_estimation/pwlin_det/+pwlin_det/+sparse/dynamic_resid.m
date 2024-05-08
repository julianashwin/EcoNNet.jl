function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = pwlin_det.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(4, 1);
    residual(1) = (y(6)) - (params(1)*y(10)+params(2)*y(5)+y(8));
    residual(2) = (y(5)) - (params(3)*y(1)-params(4)*(y(6)*params(5)-y(10))+y(7));
    residual(3) = (y(8)) - (params(8)*y(4)+x(2));
    residual(4) = (y(7)) - (params(9)*y(3)+x(1));
end
