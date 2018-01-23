clc
clear all
%% User inputs
Lx = input('length of x=');
Ly = input('length of y=');
Lz = input('length of z=');
nx = input('nx=');
ny = input('ny=');
nz = input('nz=');
ne = nx*ny*nz;
a = Lx/nx; %elemeny width
b = Ly/ny; %element length
c = Lz/nz; %element height
% Parameters
P = 0;
rho = ones(ne,1);
KE = Elementstiffness_3D2(a,b,c);
K = global_matrix3D(KE,nx,ny,nz,P,rho);
ndof = 3*(nx+1)*(ny+1)*(nz+1);
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nx, 0, 0:nz);                 % Coordinates
loadnid = kl*(nx+1)*(ny+1)+il*(ny+1)+(ny+1-jl);     % Node IDs
loaddof = 3*loadnid(:) - 1;                         % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:ny,0:nz);                  % Coordinates
fixednid = kf*(nx+1)*(ny+1)+iif*(ny+1)+(ny+1-jf);     % Node IDs
FixDOF = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
%% Force vector
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
FreeDOF = setdiff(1:ndof,FixDOF);
tolit=1e-8;
maxit=8000;
M=diag(diag(K(FreeDOF,FreeDOF)));
U(FreeDOF,:)=pcg(K(FreeDOF,FreeDOF),F(FreeDOF,:),tolit,maxit,M);
U(FixDOF,:) = 0;