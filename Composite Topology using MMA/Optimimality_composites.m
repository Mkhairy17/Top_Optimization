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
%% Parameters
P = 1;
rho_min = 10^-3;
R = 1.2*a;
eita_V0 = 0.5;
Ae = a*b;
%% Initialization of density distribution
Filter_template = ones(nx*ny,1);
iter = 1;
%% Initial A matrix values 
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
xval = 0.5*ones(5*Ne,1);
xval(Ne+1 : 2*Ne) = 0;
xval(3*Ne+1 : 4*Ne) = 0;
rho_old = xval(4*Ne+1 : 5*Ne,1);
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval,Q,Ne);
%% Non-Linear programming algorithm
figure(1)
Lamination_Parameters_Check()
figure(2)
Positivie_Definite_Check(U1,U2,U3,U4,U5)
while (1)
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval,Q,Ne);
[C,dC_drho,dC_dV1A,dC_dV2A,dC_dV3A,dC_dV4A,Strain_energy] = sensitivity_analysis_optimality_composites(rho_old,Amat,U2,U3,nx,ny,a,b,P);
f0val = C;
df0dx = [dC_dV1A;dC_dV2A;dC_dV3A;dC_dV4A;dC_drho]; 
[xmma] = Update_variables_Composites(Ne,xval,eita_V0,Ae,f0val,df0dx);
v1 = xmma(1);
v3 = xmma(3);
figure(1)
hold on
plot(v1,v3,'ro')
figure(2)
hold on
plot(v1,v3,'ro')
rho_new = xmma(4*Ne+1 : 5*Ne,1);
Error = norm(rho_old - rho_new,'inf');
%% Stopping Criteria
if  Error < 0.001
    break;
end
rho_old = rho_new;
xval = xmma;
%Lambda_old = Lambda_new;
iter = iter+1;
%% Plotting the results
x_new = reshape(rho_new,nx,ny)';
x_new = [flip(x_new')' x_new];
figure(3)
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
fprintf('Objective function %8.3f , Error in the density %8.3f , iteration %d \n',full(C),Error,iter);
end


