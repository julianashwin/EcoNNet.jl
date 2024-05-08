function [lhs, rhs] = static_resid(y, x, params)
T = NaN(0, 1);
lhs = NaN(4, 1);
rhs = NaN(4, 1);
lhs(1) = y(2);
rhs(1) = y(2)*params(1)+params(2)*y(1)+y(4);
lhs(2) = y(1);
rhs(2) = y(1)*params(3)-params(4)*(y(2)*params(5)-y(2))+y(3);
lhs(3) = y(4);
rhs(3) = y(4)*params(8)+x(2);
lhs(4) = y(3);
rhs(4) = y(3)*params(9)+x(1);
end
