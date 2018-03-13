clc
clear all
Lx = input('length of x=');
Ly = input('length of y=');
nx = input('nx=');
ny = input('ny=');
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%%
P=3;
% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = union([1:2:2*(ny+1)],[2*(nx+1)*(ny+1)]);
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%%Define forces
F(2,1)=-1;
volfrac=0.5;
rho_old=ones(ny*nx,1)*0.5;
rho_new=zeros(ny*nx,1);
rho_min=10^-3;
E=1;
v=0.3;
s_6 = [[1 -v -v;-v 1 -v;-v -v 1] zeros(3,3);zeros(3,3) 2*(1+v)*eye(3,3)]/E;
c_6 = inv(s_6);
[KE,Cm]=Elementstiffness(a,b,c_6);
%%
%Calculate initial Lambda
K=global_matrix3(KE,nx,ny,P,rho_old);
U = sparse(2*(ny+1)*(nx+1),1);
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;
strain = Calc_str(a,b,nx,ny,U);
Lambda_old=1;
iteration=1;
%%
while 1
    volfractioncalc = @(Lambda) updatedensity(Lambda,rho_old,strain,P,rho_min,Cm,a,b) /volfrac-1.0;
  options=optimset('Display','off');
    Lambda_new=fsolve(volfractioncalc,Lambda_old,options);
    [volfrac_new,rho_new]=updatedensity(Lambda_new,rho_old,strain,P,rho_min,Cm,a,b);
    norm(rho_old-rho_new,'inf')
if norm(rho_old-rho_new,'inf') < 1.0e-3
   break;
else
    rho_old=rho_new;
    Lambda_old=Lambda_new;
    K=global_matrix3(KE,nx,ny,P,rho_old);
    U = sparse(2*(ny+1)*(nx+1),1);
    U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
    U(FixDOF,:) = 0;
    strain = Calc_str(a,b,nx,ny,U);
    Objectivefunc(iteration)=full(F)'*full(U);
    iteration=iteration+1;
    plot(Objectivefunc) 
    drawnow
end
% rho_new2=reshape(rho_new,nx,ny);
%   colormap(gray); imagesc(-rho_new2'); axis equal; axis tight; axis off;pause(1e-6);
end
% figure
% rho_new2=reshape(-rho_new,nx,ny);
%   colormap(gray); imagesc(rho_new2'); axis equal; axis tight; axis off;pause(1e-6);
