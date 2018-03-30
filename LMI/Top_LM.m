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
%% Initial A matrix values 
xval_new = 0*ones(4*Ne,1);
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[Amat,U1,U2,U3,U4,U5] = xval_Amat(xval_new,Q,Ne);
%% Displacement solver 
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
[U] = Displacement_solver(Amat,nx,ny,a,b);
el = 1;
for j = 1:ny
    for i = 1:nx
       [eigen_vec] = cal_C_element(Amat,i,j,el,ny,U,KE_A11,KE_A12,KE_A16,KE_A26,KE_A22,KE_A66);
        vv1 = eigen_vec(:,1);
        vv2 = eigen_vec(:,2);
        vv3 = eigen_vec(:,3);
        xval_opt = optimum_variables (U1,U2,U3,U4,U5,vv1,vv2,vv3)
    end
end