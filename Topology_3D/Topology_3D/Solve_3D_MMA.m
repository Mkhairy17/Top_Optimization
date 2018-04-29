clc
clear all
%% User inputs
global nx ny nz Lx Ly Lz P
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
E = 1;
v = 0.3;
%% Parameters
P = 3;
rho_min = 10^-3;
Volume_Fraction_constraint = 0.2;
%% Initialization of density distribution
rho_old = ones(ny*nx*nz,1)*0.5;
iter = 1;
%% Non-Linear programming algorithm
while (1)
[C,dC] = sensitivity_analysis_3D(rho_old,nx,ny,nz,Lx,Ly,Lz,P);
rho_new = MMA_3D(nx,ny,nz,rho_old,Volume_Fraction_constraint,dC);
Error = norm(rho_old - rho_new,'inf');
%% Stopping Criteria
if  Error < 1.0e-4
    break;
end
rho_old = rho_new;
iter = iter+1;
end
