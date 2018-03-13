%---------------------------------------------------------------------
%  This is the file toymain.m.  Version September 2009.
%  Written by Krister Svanberg <krille@math.kth.se>.
%
%  This file contains a main program for using MMA to solve
%  a problem defined by the users files toyinit.m
%  (which must be run before toymain.m) and toy2.m.
%
%%%% If outeriter=0, the user should now calculate function values
%%%% and gradients of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:
if outeriter < 0.5
  [f0val,df0dx,fval,dfdx] = toy2(xval);
%  innerit=0;
  outvector1 = [outeriter xval']
  outvector2 = [f0val fval']
end
%
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
  [f0val,df0dx,fval,dfdx] = toy2(xval);
%%%% The residual vector of the KKT conditions is calculated:
  [residu,kktnorm,residumax] = ...
  kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
  outvector1 = [outeriter xval']
  outvector2 = [f0val fval']
%
end
%---------------------------------------------------------------------
