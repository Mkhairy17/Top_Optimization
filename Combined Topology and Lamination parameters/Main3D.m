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
R = 0.0001;
Lambda_old = 10^-5;
vn = Ve * ones(Ne,1);
Filter_template = ones(nx*ny*nz,1);
%=====================================
edofMat = DOF3D(nx,ny,nz);
xval_old = zeros(15*Ne,1);
rho_old = Volume_Fraction_constraint*ones(Ne,1);
xval_old(14*Ne+1 : 15*Ne) = rho_old ;
xval_new = xval_old ;
[KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66] = Matrix_derivatives(a,b,c);
iter = 1;
Compliance_old = 10^9;
 while(1)
ne = 1;
[F1111,F1122,F1133,F1123,F1113,F1112,F2222,F2233,F2223,F1223,F1222,F3333,F2333,F1333,F1233] = Lam_parameters_Fijkl(xval_old,Ne);
Cmat = Elasticity_Matrix(a1,a2,a3,a4,a5,F1111,F1122,F1133,F2222,F2233,F3333,F2333,F2223,F1123,F1113,F1233,F1333,F1112,F1222,F1223);
U = Displacement_Solver(Cmat,rho_old,a,b,c,nx,ny,nz,P);
Compliance_new = 0;
tic
parfor ne = 1 : Ne
        edof = edofMat(ne,:);
        rho_element = rho_old(ne);
%       [KE] = Elementstiffness_3D(Cmat(ne,:),a,b,c);
        [L,Phi,Compliance] = cal_C_element(Cmat,edof,U,KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66,rho_element,ne,P);
        Compliance_new = Compliance_new + Compliance;
        xval_opt = optimum_variables2(a1,a2,a3,a4,a5,P,rho_element,L);       
        V1(ne) = xval_opt(1);
        V2(ne) = xval_opt(2);
        V3(ne) = xval_opt(3);
        V4(ne) = xval_opt(4);
        V5(ne) = xval_opt(5); 
        V6(ne) = xval_opt(6); 
        V7(ne) = xval_opt(7); 
        V8(ne) = xval_opt(8); 
        V9(ne) = xval_opt(9); 
        V10(ne) = xval_opt(10);
        V11(ne) = xval_opt(11); 
        V12(ne) = xval_opt(12); 
        V13(ne) = xval_opt(13); 
        V14(ne) = xval_opt(14);

end
toc
xval_new = [V1';V2';V3';V4';V5';V6';V7';V8';V9';V10';V11';V12';V13';V14'];
if iter == 1
F = Calc_F_Sensitivities_3D(Filter_template,R,a,b,c,nx,ny,nz);
end
% Sensitivities
[rho_filtered,bn] = Filtering_subroutine_3D(U,Cmat,rho_old,F,nx,ny,nz,a,b,c,P);
[rho_new,Lambda_new] = Solve_Lambda_3D(Lambda_old,Volume_Fraction_constraint,rho_filtered,P,rho_min,nx,ny,nz,vn,bn);
rho_new = real(rho_new); 
Lambda_new = real(Lambda_new);
xval_new(14*Ne+1 : 15*Ne) = rho_new;
Error1(iter) = norm(xval_old(1:14*Ne) - xval_new(1:14*Ne),'inf');
Error2(iter) = norm(xval_old(14*Ne+1 : 15*Ne) - xval_new(14*Ne+1 : 15*Ne),'inf');
if abs(Compliance_new - Compliance_old) < 0.1 
    break;
end
Compliance_old = Compliance_new;
xval_old = xval_new;
rho_old = rho_new;
plot(iter,Compliance_new,'o')
hold on
drawnow
iter = iter + 1;
save('updated_workspace')
end
% display_3D(reshape(rho_new,ny,nx,nz))
