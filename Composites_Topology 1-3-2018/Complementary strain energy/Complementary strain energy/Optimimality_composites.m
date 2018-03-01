clc
clear all
%% User inputs
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
Ne = nx*ny;
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
vn = ones(nx*ny,1)*a*b/(Lx*Ly);
%% Parameters
P = 3;
rho_min = 10^-3;
R = 1.2*a;
eita_V0 = 0.5;
Ae = a*b;
%% Initialization of density distribution
Filter_template = ones(nx*ny,1);
iter = 1;
Lambda_old = 10;
%% Initial A matrix values 
Q = [181.81 2.90 0;2.9 10.35 0;0 0 7.71] * 1e9;
xval = 0.5*ones(5*Ne,1);
[Amat,U2,U3] = xval_Amat(xval,Q,Ne);
rho = xval(4*Ne+1 : 5*Ne,1);
%% Non-Linear programming algorithm
while (1)
%if iter == 1
%F = Calc_F_Sensitivities(Filter_template,R,a,b,nx,ny);
%end
% Filtered_densities
% [rho_old_filtered,bn,C] = Filtering_subroutine (rho_old,F,nx,ny,a,b,P);
[C,dC_drho,dC_dV1A,dC_dV2A,dC_dV3A,dC_dV4A,Strain_energy] = sensitivity_analysis_optimality_composites(rho,Amat,U2,U3,nx,ny,a,b,P);
f0val = C;
df0dx = 
[xmma] = Update_variables_Composites(Ne,xval,eita_V0,Ae,f0val,df0dx);
Error = norm(rho_old - rho_new,'inf');
%% Stopping Criteria
if  Error < 0.0001
    break;
end
rho_old = rho_new;
Lambda_old = Lambda_new;
iter = iter+1;
%% Plotting the results
x_new = reshape(rho_old_filtered,nx,ny)';
x_new = [flip(x_new')' x_new];
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
fprintf('Objective function %8.3f , Error in the density %8.3f , iteration %d \n',full(C),Error,iter);
end


