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
%===========Initializations================
P = 3;

eta_V0 = 0.5;
rho_min = 10e-3;
R = 1.4*a;
Lambda_old = 10;
vn = Ae * ones(Ne,1);
Filter_template = ones(nx*ny,1);
%=====================================
xval_old = 0.5*ones(5*Ne,1);
xval_new = xval_old ;
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
iter = 1;
while(1)
rho_old = xval_old(4*Ne+1:5*Ne);
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval_old,Q,Ne);
U_penalized = Displacement_solver_density(Amat,rho_old,nx,ny,a,b,P);
el = 1;
for j = 1:ny
    for i = 1:nx
        [eigen_vec,phi,L] = cal_C_element(Amat,i,j,el,ny,U_penalized,KE_A11,KE_A12,KE_A16,KE_A26,KE_A22,KE_A66,rho_old(el),P);
        eigen_vec = real(eigen_vec);
        vv1 = eigen_vec(:,1);
        vv2 = eigen_vec(:,2);
        vv3 = eigen_vec(:,3);
%       xval_opt = optimum_variables (U1,U2,U3,U4,U5,vv1,vv2,vv3);
%         xval_opt = optimum_variables_chol2(U1,U2,U3,U4,U5,L);
        xval_opt = optimum_variables_chol(U1,U2,U3,U4,U5,L);
        xval_new(el) = xval_opt(1);
        xval_new(el+Ne) = xval_opt(2);
        xval_new(el+2*Ne) = xval_opt(3);
        xval_new(el+3*Ne) = xval_opt(4); 
        el = el+1;
    end
end
if iter == 1
F = Calc_F_Sensitivities(Filter_template,R,a,b,nx,ny);
end
% Sensitivities
[rho_filtered,bn,C] = Filtering_subroutine (Amat,U_penalized,rho_old,F,nx,ny,a,b,P);
[rho_new,Lambda_new] = Solve_Lambda(Lambda_old,eta_V0,rho_filtered,P,rho_min,nx,ny,vn,bn);
xval_new(4*Ne+1 : 5*Ne) = rho_new;
Error1 = norm(xval_old(1:4*Ne) - xval_new(1:4*Ne),'inf')
Error2 = norm(rho_new - rho_old,'inf')
if Error1 < 0.1 && Error2 < 0.1
    break;
end
xval_old = xval_new;
iter = iter + 1;
end

%====================================================

V1 = xval_new (1:Ne);
% V2 = xval_new (Ne+1:2*Ne);
V3 = xval_new (2*Ne+1:3*Ne);
% V4 = xval_new (3*Ne+1:4*Ne);
indx = find(xval_new(4*Ne+1:5*Ne) < 1);
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
