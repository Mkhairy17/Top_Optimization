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
vn = (ones(nx*ny*nz,1)*a*b*c)/(Lx*Ly*Lz);
E = 1;
v = 0.3;
%% Parameters
P = 3;
rho_min = 10^-3;
R = 2*a;
Volume_Fraction_constraint = 0.5;
%% Initialization of density distribution
rho_old = ones(ny*nx*nz,1)*0.5;
Filter_template = ones(ny*nx*nz,1);
iter = 1;
Lambda_old = 10;
%% Non-Linear programming algorithm
while (1)
if iter == 1
F = Calc_F_Sensitivities_3D(Filter_template,R,a,b,c,nx,ny,nz);
end
% Filtered_densities
[rho_old_filtered,bn,C] = Filtering_subroutine_3D(rho_old,F,nx,ny,nz,Lx,Ly,Lz,P);
[rho_new,Lambda_new] = Solve_Lambda_3D(Lambda_old,Volume_Fraction_constraint,rho_old,P,rho_min,nx,ny,nz,vn,bn);
Error = norm(rho_old - rho_new,'inf')
%% Stopping Criteria
if  Error < 0.01
    break;
end
rho_old = rho_new;
iter = iter+1;
end
%%
display_3D(reshape(rho_old_filtered,ny,nx,nz))