function [lhs, rhs] = static_resid(y, x, params)
T = NaN(0, 1);
lhs = NaN(5, 1);
rhs = NaN(5, 1);
lhs(1) = y(2);
rhs(1) = params(1)*y(3)+params(2)*y(1)+y(5);
lhs(2) = y(2)-y(3);
rhs(2) = x(3);
lhs(3) = y(1);
rhs(3) = y(1)*params(3)-params(4)*(y(2)*params(5)-y(3))+y(4);
lhs(4) = y(5);
rhs(4) = y(5)*params(6)+x(2);
lhs(5) = y(4);
rhs(5) = y(4)*params(7)+x(1);
end
