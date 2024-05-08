function [lhs, rhs] = dynamic_resid(y, x, params, steady_state)
T = NaN(0, 1);
lhs = NaN(5, 1);
rhs = NaN(5, 1);
lhs(1) = y(7);
rhs(1) = params(1)*y(8)+params(2)*y(6)+y(10);
lhs(2) = y(7)-y(3);
rhs(2) = x(3);
lhs(3) = y(6);
rhs(3) = params(3)*y(1)-params(4)*(y(7)*params(5)-y(8))+y(9);
lhs(4) = y(10);
rhs(4) = params(6)*y(5)+x(2);
lhs(5) = y(9);
rhs(5) = params(7)*y(4)+x(1);
end
