%---------------------------------------------------------------------
%  This is the file toy2.m.  Version September 2009.
%  Written by Krister Svanberg <krille@math.kth.se>.
%  It calculates function values and gradients
%  for the following "toy problem":
%
%    minimize x(1)^2 + x(2)^2 + x(3)^2
%  subject to (x(1)-5)^2 + (x(2)-2)^2 + (x(3)-1)^2 =< 9
%             (x(1)-3)^2 + (x(2)-4)^2 + (x(3)-3)^2 =< 9
%              0 =< x(j) =< 5, for j=1,2,3.
%
%  (the bounds on x(j) are defined in toyinit.m)
%
function [f0val,df0dx,fval,dfdx] = toy2(x);
%
f0val = x(1)^2 + x(2)^2 + x(3)^2;
%
df0dx = [2*x(1)
	 2*x(2)
	 2*x(3)];
%
fval  = [(x(1)-5)^2+(x(2)-2)^2+(x(3)-1)^2-9
	 (x(1)-3)^2+(x(2)-4)^2+(x(3)-3)^2-9];
%
dfdx  = 2*[x(1)-5  x(2)-2  x(3)-1
	   x(1)-3  x(2)-4  x(3)-3];
%---------------------------------------------------------------------
