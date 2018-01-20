clc
clear all
%% User inputs
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
ne = nx*ny;
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%% Parameters
P = 3;
rho_min = 10^-3;
Volume_Fraction_constraint = 0.3;
R = 0.04;
E = 1;
v = 0.3;
s_6 = [[1 -v -v;-v 1 -v;-v -v 1] zeros(3,3);zeros(3,3) 2*(1+v)*eye(3,3)]/E;
c_6 = inv(s_6);
[KE,Cm] = Elementstiffness(a,b,c_6);
% Initialization of density distribution
rho_old = ones(ny*nx,1)*0.5;
iter = 1;
[strain] = Find_Strain(KE,rho_old,Lx,Ly,nx,ny,P);
 %%
%Calculate initial Lambda
[Strain_energy_new,Lambda0] = Calc_Strain_Energy(rho_old,strain,R,P,nx,ny,a,b,Cm,0);
%% Optimiality criteria
while (1)
strain = Find_Strain(KE,rho_old,Lx,Ly,nx,ny,P);
if iter == 1
F = Calc_F(rho_old,strain,R,a,b,nx,ny,P,Cm);
end
[rho_new,Lambda_new] = Solve_Volume_Fraction(Lambda0,Volume_Fraction_constraint,rho_old,P,rho_min,R,a,b,nx,ny,strain,Cm,F);
[Strain_energy_new] = Calc_Strain_Energy(rho_old,strain,R,P,nx,ny,a,b,Cm,F);
Total_strain_energy = sum(Strain_energy_new);
Error = norm(rho_old-rho_new,'inf');
%% Stopping Criteria
if  Error < 1.0e-5
    break;
end
rho_old = rho_new;
Lambda0 = Lambda_new;

%% Plotting the results
x_new = reshape(rho_new,nx,ny)';
colormap(gray); imagesc(-x_new); axis equal; axis tight; axis off;pause(1e-6);
iter = iter+1;
fprintf('Objective function %8.3f , Error in the density %8.3f , Lambda %8.3f iteration %d \n',Total_strain_energy,Error,Lambda_new,iter);
end


