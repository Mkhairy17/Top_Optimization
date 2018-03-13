%---------------------------------------------------------------------
%  This is the file trussinit.m
%  in which some vectors for the "three bar truss problem",
%  defined in the file trussmain.m, are initialized.
%
%  Written by Krister Svanberg <krille@math.kth.se>
%  Department of Mathematics
%  SE-10044 Stockholm, Sweden.
%
m = 4;
n = 3;
epsimin = 0.0000001;
xval  = ones(n,1);
xold1 = xval;
xold2 = xval;
xmin  = 0.001*ones(n,1);
xmax  = 3*ones(n,1);
low   = xmin;
upp   = xmax;
c = 1000*ones(m,1);
d = zeros(m,1);
a0 = 1;
a = [1 1 1 0]';
%
outeriter = 0;
maxoutit  = 1;
kkttol  = 0;
%
%---------------------------------------------------------------------

