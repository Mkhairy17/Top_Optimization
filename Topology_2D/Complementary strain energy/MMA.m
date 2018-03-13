function [xnew]=MMA(nelx,nely,x,volfrac,dc)
xlow=0.001; xhigh=1;
L=x-0.1*(xhigh-xlow)*ones(nely,nelx);
high=(x-L).^2.*-dc./(xlow-L).^2;
low=(x-L).^2.*-dc./(xhigh-L).^2; 
l2 = min(min(high));
l1 = max(max(low));
for i=1:50
lmid = 0.5*(l2+l1);
xnew=max(xlow,min(xhigh,L+abs(x-L).*sqrt(-dc./lmid)));
if sum(sum(xnew)) - volfrac*nelx*nely > 0;
l2 = lmid;
else
l1 = lmid;
end
end