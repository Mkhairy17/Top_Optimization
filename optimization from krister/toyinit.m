%---------------------------------------------------------------------
%  This is the file toyinit.m.  Version September 2009.
%  Written by Krister Svanberg <krille@math.kth.se>
%  Some parameters and the starting point are defined
%  for the "toy problem":
%    minimize x(1)^2 + x(2)^2 + x(3)^2
%  subject to (x(1)-5)^2 + (x(2)-2)^2 + (x(3)-1)^2 =< 9
%             (x(1)-3)^2 + (x(2)-4)^2 + (x(3)-3)^2 =< 9
%              0 =< x(j) =< 5, for j=1,2,3.
%
%
m = 2;
n = 3;
epsimin = 0.0000001;
xval    = [4 3 2]';
xold1   = xval;
xold2   = xval;
xmin    = [0  0  0]';
xmax    = [5  5  5]';
low     = xmin;
upp     = xmax;
c       = [1000  1000]';
d       = [1  1]';
a0      = 1;
a       = [0  0]';;
outeriter = 0;
maxoutit  = 1;
kkttol  = 0;
%
%---------------------------------------------------------------------
