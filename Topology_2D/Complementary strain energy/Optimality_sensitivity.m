clc
clear all
%% User inputs
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
ne = nx*ny;
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
vn = ones(nx*ny,1)*a*b;
%% Parameters
P = 3;
rho_min = 10^-3;
Volume_Fraction_constraint = 0.5;
%% Initialization of density distribution
rho_old = ones(ny*nx,1)*0.5;
iter = 1;
Lambda_old=10;
%% Non-Linear programming algorithm
while (1)
[C,dC,Strain_energy,bn] = sensitivity_analysis_optimality(rho_old,nx,ny,Lx,Ly,P);
[rho_new,Lambda_new] =Solve_Lambda(Lambda_old,Volume_Fraction_constraint,rho_old,P,rho_min,nx,ny,vn,bn);
Error = norm(rho_old - rho_new,'inf');
%% Stopping Criteria
if  Error < 0.01
    break;
end
rho_old = rho_new;
Lambda_old=Lambda_new;
iter = iter+1;
%% Plotting the results
x_new = reshape(rho_new,nx,ny)';
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
fprintf('Objective function %8.3f , Error in the density %8.3f , iteration %d \n',full(C),Error,iter);
end


