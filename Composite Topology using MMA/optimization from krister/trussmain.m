%---------------------------------------------------------------------
%  This is the file trussmain.m.  Version September 2009.
%  Written by Krister Svanberg <krille@math.kth.se>.
%
%  This file contains a main program for using MMA to solve
%  a problem defined by the users files trussinit.m
%  (which must be run before trussmain.m) and truss2.m.
%
%  The following "three bar truss problem" is solved:
%
%  mininimize max{f1(x),f2(x),f3(x)}
%  subject to:
%  x1 + x2 + x3 <= V,
%   0.001 <= xj <= V, for j=1,2,3.
%
%  where
%  xj = volume of the j:th bar,
%  fi(x) = pi'*ui = compliance for the i:th loadcase,
%  V = 3 = upper bound on the total volume.
%
%  Since the length of each bar is = 1,
%  xj is also the cross section area of the j:the bar.
%  
%  Bar 1 connects the nodes 1 and 4.
%  Bar 2 connects the nodes 2 and 4.
%  Bar 3 connects the nodes 3 and 4.
%  The coordinates of node 1 are (-1,0).
%  The coordinates of node 2 are (-1/sqrt(2),-1/sqrt(2)).
%  The coordinates of node 3 are (0,-1).
%  The coordinates of node 4 are (0,0).
%  The nodes 1,2,3 are fixed, while the load vectors for
%  the different loadcases are applied at node 4.
%  The load vector p1 = (1,0)' (loadcase 1).
%  The load vector p2 = (1,1)' (loadcase 2).
%  The load vector p3 = (0,1)' (loadcase 3).
%
%  The displacement vector ui is obtained from the system
%  K(x)*ui = pi, i=1,2,3, where K(x) is the stiffness matrix
%  and pi is the load vector.
%  The stiffness matrix is given by K(x) = R*D(x)*R',
%  where D(x) is a diagonal matrix with diagonal elements x1,x2,x3.
%  The constant matrix R is given in the code below.
%  The derivatives of the functions fi(x) are then given by
%  dfi/dxj = -(rj'*ui)^2 ,
%  where rj is the j:th column of the matrix R.
%
%  The problem is written on the following form required by MMA.
%
%  minimize  z + 1000*(y1+y2+y3+y4)
%  subject to the constraints:
%            f1(x) - z - y1 <= 0
%            f2(x) - z - y2 <= 0
%            f3(x) - z - y3 <= 0
%     x1 + x2 + x3 - 3 - y4 <= 0
%                   0 <= xj <= 3, j=1,2,3
%                        yi >= 0, i=1,2,3,4
%                         z >= 0.
%
%
%%%% If outeriter=0, the user should now calculate function values
%%%% and gradients of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:
if outeriter < 0.5
  [f0val,df0dx,fval,dfdx] = truss2(xval);
  innerit=0;
  outvector1 = [outeriter xval']
  outvector2 = fval'
end
%%%% The iterations start:
kktnorm = kkttol+10;
outit = 0;
while kktnorm > kkttol & outit < maxoutit
  outit   = outit+1;
  outeriter = outeriter+1;
%%%% The MMA subproblem is solved at the point xval:
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
%%%% Some vectors are updated:
  xold2 = xold1;
  xold1 = xval;
  xval  = xmma;
%%%% The user should now calculate function values and gradients
%%%% of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx.
  [f0val,df0dx,fval,dfdx] = truss2(xval);
%%%% The residual vector of the KKT conditions is calculated:
  [residu,kktnorm,residumax] = ...
  kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
  outvector1 = [outeriter xval']
  outvector2 = fval'
%
end
%---------------------------------------------------------------------
