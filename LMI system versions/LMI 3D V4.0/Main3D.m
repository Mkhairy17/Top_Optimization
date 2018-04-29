clc
clear all
%% User inputs
Lx = input('length of x=');
Ly = input('length of y=');
Lz = input('length of z=');
nx = input('nx=');
ny = input('ny=');
nz = input('nz=');
Ne = nx*ny*nz;
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
c = Lz/nz;
Ve = a*b*c;
%===========Initializations================
a1 = 1;
a2 = 1;
a3 = 1;
a4 = 1;
a5 = 1;
P = 3;
Volume_Fraction_constraint = 0.5;
rho_min = 10e-3;
R = a;
Lambda_old = 10;
vn = Ve * ones(Ne,1);
Filter_template = ones(nx*ny*nz,1);
%=====================================
edofMat = DOF3D(nx,ny,nz);
xval_old = 0*ones(14*Ne,1);
rho_old = ones(Ne,1);
xval_old(14*Ne+1 : 15*Ne) = rho_old ;
xval_new = xval_old ;
[KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66] = Matrix_derivatives(a,b,c);
iter = 1;
while(1)
ne = 1;
[F1111,F1122,F1133,F1123,F1113,F1112,F2222,F2233,F2223,F1223,F1222,F3333,F2333,F1333,F1233] = Lam_parameters_Fijkl(xval_old,Ne);
Cmat = Elasticity_Matrix(a1,a2,a3,a4,a5,F1111,F1122,F1133,F2222,F2233,F3333,F2333,F2223,F1123,F1113,F1233,F1333,F1112,F1222,F1223);
U = Displacement_Solver(Cmat,rho_old,a,b,c,nx,ny,nz,P);
for j = 1:ny
    for i = 1:nx
        for k = 1:nz
        edof = edofMat(ne,:);
        rho_element = rho_old(ne);
        [L,Phi] = cal_C_element(Cmat,edof,U,KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66,rho_element,ne,P);
        xval_opt = optimum_variables2(a1,a2,a3,a4,a5,L);       
        xval_new(ne) = xval_opt(1);
        xval_new(ne+Ne) = xval_opt(2);
        xval_new(ne+2*Ne) = xval_opt(3);
        xval_new(ne+3*Ne) = xval_opt(4);
        xval_new(ne+4*Ne) = xval_opt(5); 
        xval_new(ne+5*Ne) = xval_opt(6); 
        xval_new(ne+6*Ne) = xval_opt(7); 
        xval_new(ne+7*Ne) = xval_opt(8); 
        xval_new(ne+8*Ne) = xval_opt(9); 
        xval_new(ne+9*Ne) = xval_opt(10);
        xval_new(ne+10*Ne) = xval_opt(11); 
        xval_new(ne+11*Ne) = xval_opt(12); 
        xval_new(ne+12*Ne) = xval_opt(13); 
        xval_new(ne+13*Ne) = xval_opt(14); 
        ne = ne+1;
        end
    end
end
if iter == 1
F = Calc_F_Sensitivities_3D(Filter_template,R,a,b,c,nx,ny,nz);
end
% Sensitivities
[rho_filtered,bn] = Filtering_subroutine_3D(U,Cmat,rho_old,F,nx,ny,nz,a,b,c,P);
[rho_new,Lambda_new] = Solve_Lambda_3D(Lambda_old,Volume_Fraction_constraint,rho_old,P,rho_min,nx,ny,nz,vn,bn);
xval_new(14*Ne+1 : 15*Ne) = rho_new;
Error = norm(xval_old - xval_new,'inf')
if Error < 0.1 
    break;
end
xval_old = xval_new;
rho_old = rho_new;
iter = iter + 1;
end

% %====================================================
% 
% V1 = xval_new (1:Ne);
% % V2 = xval_new (Ne+1:2*Ne);
% V3 = xval_new (2*Ne+1:3*Ne);
% % V4 = xval_new (3*Ne+1:4*Ne);
% indx = find(xval_new(4*Ne+1:5*Ne) < 1);
% figure(1)
% V1(indx) = -1.1;
% V1_reshaped = reshape(V1,nx,ny)';
% imagesc(V1_reshaped)
% colorbar
% colormap jet;
% title('V1')
% figure(2)
% V3(indx) = -1.1;
% V3_reshaped = reshape(V3,nx,ny)';
% imagesc(V3_reshaped)
% colorbar
% colormap jet;
% title('V3')
% 
% 
% % V2_reshaped = reshape(V2,nx,ny)';
% % V4_reshaped = reshape(V4,nx,ny)';
