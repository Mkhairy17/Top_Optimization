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
eta_V0 = 0.5;
rho_min = 10e-3;
R = 1.4*a;
Lambda_old = 10;
vn = Ve * ones(Ne,1);
Filter_template = ones(nx*ny,1);
%=====================================
edofMat = DOF3D(nx,ny,nz);
xval_old = 0*ones(14*Ne,1);
xval_new = xval_old ;
[KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66] = Matrix_derivatives(a,b,c);
iter = 1;
while(1)
ne = 1;
rho_old = ones(Ne,1);
[F1111,F1122,F1133,F1123,F1113,F1112,F2222,F2233,F2223,F1223,F1222,F3333,F2333,F1333,F1233] = Lam_parameters_Fijkl(xval_old,Ne);
Cmat = Elasticity_Matrix(a1,a2,a3,a4,a5,F1111,F1122,F1133,F2222,F2233,F3333,F2333,F2223,F1123,F1113,F1233,F1333,F1112,F1222,F1223);
U = Displacement_Solver(Cmat,rho_old,a,b,c,nx,ny,nz,P);
for j = 1:ny
    for i = 1:nx
        for k = 1:nz
        edof = edofMat(ne,:);
        rho_element = rho_old(ne);
        [L,Phi] = cal_C_element(Cmat,edof,U,KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66,rho_element,ne,P);
        xval_opt = optimum_variables (a1,a2,a3,a4,a5,L);       
        xval_new(ne) = xval_opt(22);
        xval_new(ne+Ne) = xval_opt(23);
        xval_new(ne+2*Ne) = xval_opt(24);
        xval_new(ne+3*Ne) = xval_opt(25);
        xval_new(ne+4*Ne) = xval_opt(26); 
        xval_new(ne+5*Ne) = xval_opt(27); 
        xval_new(ne+6*Ne) = xval_opt(28); 
        xval_new(ne+7*Ne) = xval_opt(29); 
        xval_new(ne+8*Ne) = xval_opt(30); 
        xval_new(ne+9*Ne) = xval_opt(31);
        xval_new(ne+10*Ne) = xval_opt(32); 
        xval_new(ne+11*Ne) = xval_opt(33); 
        xval_new(ne+12*Ne) = xval_opt(34); 
        xval_new(ne+13*Ne) = xval_opt(35); 
        ne = ne+1;
        end
    end
end
Error1 = norm(xval_old(1:14*Ne) - xval_new(1:14*Ne),'inf')
if Error1 < 0.1 
    break;
end
xval_old = xval_new;
iter = iter + 1;
end


% if iter == 1
% F = Calc_F_Sensitivities(Filter_template,R,a,b,nx,ny);
% end
% % Sensitivities
% [rho_filtered,bn,C] = Filtering_subroutine (Amat,U_penalized,rho_old,F,nx,ny,a,b,P);
% [rho_new,Lambda_new] = Solve_Lambda(Lambda_old,eta_V0,rho_filtered,P,rho_min,nx,ny,vn,bn);
% xval_new(4*Ne+1 : 5*Ne) = rho_new;
% Error1 = norm(xval_old(1:4*Ne) - xval_new(1:4*Ne),'inf')
% Error2 = norm(rho_new - rho_old,'inf')
% if Error1 < 0.1 && Error2 < 0.1
%     break;
% end
% xval_old = xval_new;
% iter = iter + 1;
% end
% 
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
