function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = illus_det_est_learn_martin.static_g1_tt(T, y, x, params);
end
g1 = zeros(4, 4);
g1(1,1)=(-params(2));
g1(1,2)=1-params(1);
g1(1,4)=(-1);
g1(2,1)=1-params(3);
g1(2,2)=params(4)*(params(5)-1);
g1(2,3)=(-1);
g1(3,4)=1-params(8);
g1(4,3)=1-params(9);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
