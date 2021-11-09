function [nzij_pred, nzij_current, nzij_fwrd] = dynamic_g1_nz()
% Returns the coordinates of non-zero elements in the Jacobian, in column-major order, for each lead/lag (only for endogenous)
  nzij_pred = zeros(2, 2, 'int32');
  nzij_pred(1,1)=3; nzij_pred(1,2)=1;
  nzij_pred(2,1)=2; nzij_pred(2,2)=3;
  nzij_current = zeros(6, 2, 'int32');
  nzij_current(1,1)=1; nzij_current(1,2)=1;
  nzij_current(2,1)=3; nzij_current(2,2)=1;
  nzij_current(3,1)=1; nzij_current(3,2)=2;
  nzij_current(4,1)=2; nzij_current(4,2)=2;
  nzij_current(5,1)=3; nzij_current(5,2)=2;
  nzij_current(6,1)=1; nzij_current(6,2)=3;
  nzij_fwrd = zeros(0, 2, 'int32');
end
