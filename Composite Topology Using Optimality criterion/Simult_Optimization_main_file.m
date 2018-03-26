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
Lambda_old = 10;
%% Initial A matrix values 
xval_old = 0*ones(4*Ne,1);
xval_old(4*Ne+1:5*Ne) = 1;
xval_new = zeros(5*Ne,1);
rho_old = ones(Ne,1);
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
%% Displacement solver 
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
indx = 1;
%% Lamination parameters optimization
while(1)
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval_old,Q,Ne,P);
[U] = Displacement_solver(Amat,nx,ny,a,b);
el = 1;
for j = 1:ny
    for i = 1:nx
    xval = xval_old([el,el+Ne,el+2*Ne,el+3*Ne,el+4*Ne]);    
    cal_C_element(xval);  
    x = Varivbles_Update_Composite(xval);
    xval_new(el) = x(1);
    xval_new(el+Ne) = x(2);
    xval_new(el+2*Ne) = x(3);
    xval_new(el+3*Ne) = x(4); 
    xval_new(el+4*Ne) = x(5);   
    el = el + 1;
    end
end
Error = norm(xval_new - xval_old,'inf')
if Error < 0.01
    break;
end
xval_old = xval_new;
end

%% Material Lay-out Optimization
while (1)
[Amat_new,U1,U2,U3,U4,U5] = xval_Amat(xval_new,Q,Ne,P);
% Displacement solver 
U_new = Displacement_solver(Amat_new,nx,ny,a,b);
% Sensitivities
bn = bn_calculation (U_new,Amat_new,rho_old,a,b,nx,ny,P);
[rho_new,Lambda_new] = Solve_Lambda(Lambda_old,eta_V0,rho_old,P,rho_min,nx,ny,vn,bn);
xval_new(4*Ne+1 : 5*Ne) = rho_new;
eta = sum(rho_new) * Ae;
Error = abs(eta - eta_V0);
if Error < 0.01
    break;
end
Lambda_old = Lambda_new;
x_new = reshape(rho_new,nx,ny)';
x_new = [flip(x_new')' x_new];
figure(3)
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
end