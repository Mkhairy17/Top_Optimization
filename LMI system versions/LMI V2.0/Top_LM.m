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
Ae = a*b;

%====================================================
P = 3;
eta_V0 = 0.5;
rho_min = 10e-3;
R = 3*a;
Lambda_old = 10;
rho_old = ones(Ne,1);
vn = Ae * ones(Ne,1);
Filter_template = ones(nx*ny,1);
iter = 1;
%% Initial A matrix values 
xval_old = 0*ones(4*Ne,1);
xval_new = xval_old ;
xval_new(4*Ne+1 : 5*Ne) = rho_old;
alfa = 1*ones(3*Ne,1);
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
%% Displacement solver 
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
while(1)
    if iter == 1
F = Calc_F_Sensitivities(Filter_template,R,a,b,nx,ny);
    end
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval_old,Q,Ne);
U = Displacement_solver(Amat,[],nx,ny,a,b); %nonpenalized stiffness matrix
el = 1;
for j = 1:ny
    for i = 1:nx
        [eigen_vec,phi] = cal_C_element(Amat,i,j,el,ny,U,KE_A11,KE_A12,KE_A16,KE_A26,KE_A22,KE_A66);
        vv1 = eigen_vec(:,1);
        vv2 = eigen_vec(:,2);
        vv3 = eigen_vec(:,3);
        xval_opt = optimum_variables (U1,U2,U3,U4,U5,vv1,vv2,vv3);
        alfa(el) = xval_opt(1);
        alfa(el+1) = xval_opt(2);
        alfa(el+2) = xval_opt(3);
        xval_new(el) = xval_opt(4);
        xval_new(el+Ne) = xval_opt(5);
        xval_new(el+2*Ne) = xval_opt(6);
        xval_new(el+3*Ne) = xval_opt(7); 
        el = el+1;
    end
end
Error = norm(xval_new - xval_old,'inf')
if Error < 0.01
    break;
end
xval_old = xval_new;
end

while (1)
% Material Lay-out Optimization
Amat_new_density = xval_Amat(xval_new,Q,Ne);
% Displacement solver 
U_penalized = Displacement_solver_density(Amat_new_density,rho_old,nx,ny,a,b,P); %Penalized stiffness values
% Sensitivities
% bn = bn_calculation (U_penalized,Amat_new_density,rho_old,a,b,nx,ny,P);
[rho_filtered,bn,C] = Filtering_subroutine (Amat_new_density,U_penalized,rho_old,F,nx,ny,a,b,P);
[rho_new,Lambda_new] = Solve_Lambda(Lambda_old,eta_V0,rho_old,P,rho_min,nx,ny,vn,bn);
xval_new(4*Ne+1 : 5*Ne) = rho_new;
Error = norm(rho_new - rho_old,'inf');
if Error < 0.01
    break;
end
rho_old = rho_new;
xval_new(4*Ne+1 : 5*Ne) = rho_old;
x_new = reshape(rho_new,nx,ny)';
% x_new = [flip(x_new')' x_new];
figure(3)
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
iter = iter+1;
end

V1 = xval_new (1:Ne);
% V2 = xval_new (Ne+1:2*Ne);
V3 = xval_new (2*Ne+1:3*Ne);
% V4 = xval_new (3*Ne+1:4*Ne);
indx = find(rho_new < 1);
figure(1)
V1(indx) = -1.1;
V1_reshaped = reshape(V1,nx,ny)';
imagesc(V1_reshaped)
colorbar
colormap jet;
title('V1')
figure(2)
V3(indx) = -1.1;
V3_reshaped = reshape(V3,nx,ny)';
imagesc(V3_reshaped)
colorbar
colormap jet;
title('V3')


% V2_reshaped = reshape(V2,nx,ny)';
% V4_reshaped = reshape(V4,nx,ny)';
