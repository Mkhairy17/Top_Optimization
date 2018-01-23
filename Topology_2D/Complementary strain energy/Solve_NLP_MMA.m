clc
clear all
%% User inputs
global nx ny Lx Ly P
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
ne = nx*ny;
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%% Parameters
P = 3;
rho_min = 10^-3;
Volume_Fraction_constraint = 0.5;
%% Initialization of density distribution
rho_old = ones(ny*nx,1)*0.5;
iter = 1;

%% Non-Linear programming algorithm
while (1)
[C,dC] = sensitivity_analysis(rho_old);
rho_old = reshape(rho_old,nx,ny)';
dC = reshape(dC,nx,ny)';
rho_new = MMA(nx,ny,rho_old,0.5,dC);
Total_strain_energy = sum(full(C));
Error = norm(rho_old - rho_new,'inf');
%% Stopping Criteria
if  Error < 1.0e-4
    break;
end
rho_new = reshape(rho_new,nx*ny,1);
rho_old = rho_new;
iter = iter+1;
%% Plotting the results
x_new = reshape(rho_new,nx,ny)';
colormap(gray); imagesc(-x_new'); axis equal; axis tight; axis off;pause(1e-6);
fprintf('Objective function %8.3f , Error in the density %8.3f , iteration %d \n',Total_strain_energy,Error,iter);
end


