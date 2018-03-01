function [xmma] = Update_variables_Composites(Ne,xval,eita_V0,Ae,f0val,df0dx)
n = 5*Ne;
m = 2*Ne+1;
[dfdx]= dfdx_calc(xval,Ne,Ae);
[fval] = fval_MMA(xval,eita_V0,Ae,Ne);
xold1 = xval;
xold2 = xval;
xmin = [];
xmax = [];
low = [];
upp = [];
iter = 1;
[xmma] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,[],[],[],[]);