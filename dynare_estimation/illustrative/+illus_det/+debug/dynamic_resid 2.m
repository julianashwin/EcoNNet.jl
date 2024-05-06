function [lhs, rhs] = dynamic_resid(y, x, params, steady_state)
T = NaN(0, 1);
lhs = NaN(4, 1);
rhs = NaN(4, 1);
lhs(1) = y(6);
rhs(1) = params(1)*y(10)+params(2)*y(5)+y(8);
lhs(2) = y(5);
rhs(2) = params(3)*y(1)-params(4)*(y(6)*params(5)-y(10))+y(7);
lhs(3) = y(8);
rhs(3) = params(8)*y(4)+x(2);
lhs(4) = y(7);
rhs(4) = params(9)*y(3)+x(1);
end
