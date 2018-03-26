clc
clear all
%% User inputs
global Ne nx ny a b i j el P Amat U  Ae  Lambda KE_A11 KE_A22 KE_A66 KE_A12 KE_A16 KE_A26
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
eta_V0 = 0.4;
Ae = a*b;
vn = Ae*ones(Ne,1);
%% Initialization of Optimization multiplier
Lambda = 10;
%% Initial A matrix values 
xval_new = 0*ones(5*Ne,1);
xval_new(4*Ne+1:5*Ne) = 1;
rho_old = ones(Ne,1);
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval_new,Q,Ne,P);
%% Displacement solver 
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
[U] = Displacement_solver(Amat,nx,ny,a,b);
indx = 1;
%% Lamination parameters optimization
while(1)
el = 1;
for j = 1:ny
    for i = 1:nx
    xval = xval_new([el,el+Ne,el+2*Ne,el+3*Ne,el+4*Ne]);  
    cal_C_element(xval);  
    x = Varivbles_Update_Composit e_old(xval);
    xval_new(el) = x(1);
    xval_new(el+Ne) = x(2);
    xval_new(el+2*Ne) = x(3);
    xval_new(el+3*Ne) = x(4); 
    xval_new(el+4*Ne) = x(5);   
    el = el + 1;
    end
end
rho_new = xval_new(4*Ne+1 : 5*Ne);
eta = sum(rho_new) * Ae
Lambda = Lambda *(1 + (P+1) * (1- eta_V0/eta))
%% Plotting the results
% x_new = reshape(rho_new,nx,ny)';
%x_new = [flip(x_new')' x_new];
% figure(3)
% colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
% fprintf('Objective function %8.3f , Error in the density %8.3f , iteration %d \n',full(C),Error,iter);
% break;
end