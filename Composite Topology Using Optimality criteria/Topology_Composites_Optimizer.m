clc
clear all
%% User inputs
global Ne nx ny a b i j el P rho Amat U  Ae  Lambda
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
Ne = nx*ny;
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%% Parameters
P = 1;
rho_min = 10^-3;
eita_V0 = 0.5;
Ae = a*b;
%% Initialization of density distribution
iter = 1;
Lambda = 1000;
%% Initial A matrix values 
xval = 0.5*ones(5,1);
rho = 0.5*ones(Ne,1);
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval,Q,Ne);
%% Non-Linear programming algorithm
while (1)
[U] = Displacement_solver(rho,Amat,nx,ny,a,b,P);
el = 1;
for j = 1:ny
for i = 1:nx
    [C_element,Psi] = cal_C_element(xval);
    x = Varivbles_Update_Composite(xval);
    el = el + 1;
end
end
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval,Q,Ne);
[C,dC_drho,dC_dV1A,dC_dV2A,dC_dV3A,dC_dV4A,Strain_energy] = sensitivity_analysis_optimality_composites(rho_old,Amat,U2,U3,nx,ny,a,b,P);
f0val = C;
df0dx = [dC_dV1A;dC_dV2A;dC_dV3A;dC_dV4A;dC_drho]; 
[xmma] = Update_variables_Composites(Ne,xval,eita_V0,Ae,f0val,df0dx);
rho_new = xmma(4*Ne+1 : 5*Ne,1);
Error = norm(rho_old - rho_new,'inf')
%% Stopping Criteria
% if  Error < 0.001
%     break;
% end
rho_old = rho_new;
xval = xmma;
%Lambda_old = Lambda_new;
iter = iter+1;
%% Plotting the results
% x_new = reshape(rho_new,nx,ny)';
% x_new = [flip(x_new')' x_new];
% figure(3)
% colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
% fprintf('Objective function %8.3f , Error in the density %8.3f , iteration %d \n',full(C),Error,iter);
end


