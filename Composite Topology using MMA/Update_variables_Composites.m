function [xmma] = Update_variables_Composites(Ne,xval,eita_V0,Ae,f0val,df0dx)
n = 5*Ne;
m = 2*Ne+1;
[dfdx]= dfdx_calc(xval,Ne,Ae);
[fval] = fval_MMA(xval,eita_V0,Ae,Ne);
xold1 = xval;
xold2 = xval;
xmin = -1*ones(5*Ne,1);
xmax = 1*ones(5*Ne,1);
xmax (4*Ne+1 : 5*Ne,1) = 1;

% xmin(Ne+1:2*Ne,1) = 0;
% xmax(Ne+1:2*Ne,1) = 0.1;
% xmin(3*Ne+1:4*Ne,1) = 0;
% xmax(3*Ne+1:4*Ne,1) = 0.1;

low = -1*ones(5*Ne,1);
upp = 1*ones(5*Ne,1);
upp (4*Ne+1 : 5*Ne,1) = 1;
% low(Ne+1:2*Ne,1) = 0;
% upp(Ne+1:2*Ne,1) = 0.1;
% low(3*Ne+1:4*Ne,1) = 0;
% upp(3*Ne+1:4*Ne,1) = 0.1;

iter = 1;
[xmma] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,0,zeros(m,1),zeros(m,1),zeros(m,1));