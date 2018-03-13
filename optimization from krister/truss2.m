%---------------------------------------------------------------------
%  This is the file truss2.m
%  which calculates function values and derivatives for the
%  "three bar truss problem", defined in the file trussmain.m.
%  Written by Krister Svanberg <krille@math.kth.se>
%  Department of Mathematics
%  SE-10044 Stockholm, Sweden.
%
function [f0val,df0dx,fval,dfdx] = truss2(x);
%
e = [1 1 1]';
f0val = 0;
df0dx = 0*e;

D = diag(x);
sq2 = 1/sqrt(2);
R = [ 1 sq2 0
      0 sq2 1 ];

p1 =  [1 0]';
p2 =  [1 1]';
p3 =  [0 1]';

K = R*D*R';
u1 = K\p1;
u2 = K\p2;
u3 = K\p3;

compl1 = p1'*u1;
compl2 = p2'*u2;
compl3 = p3'*u3;

volume = e'*x;
V = 3;
vol1 = volume - V;

fval = [compl1 compl2 compl3 vol1]';

rtu1 = R'*u1;
rtu2 = R'*u2;
rtu3 = R'*u3;

dcompl1 = -rtu1.*rtu1;
dcompl2 = -rtu2.*rtu2;
dcompl3 = -rtu3.*rtu3;

dfdx = [dcompl1 dcompl2 dcompl3 e]';

%---------------------------------------------------------------------
