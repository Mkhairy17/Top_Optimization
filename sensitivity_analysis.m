function [C,dC,Volume_fraction, dVolume_fraction] = sensitivity_analysis(Lx,Ly,nx,ny,P)
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
%% Parameters
E = 1;
v = 0.3;
s_6 = [[1 -v -v;-v 1 -v;-v -v 1] zeros(3,3);zeros(3,3) 2*(1+v)*eye(3,3)]/E;
c_6 = inv(s_6);
[KE,Cm] = Elementstiffness(a,b,c_6);
% Initialization of density distribution
rho_old = ones(ny*nx,1)*0.5;
 %%
% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = union([1:2:2*(ny+1)],[2*(nx+1)*(ny+1)]);
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%%Define forces
F(2,1) = -1;
%% Global stiffness matrix
K = global_matrix3 (KE,nx,ny,P,rho_old);
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;
ne=1;
   for j = 1:ny
   for i = 1:nx
       n1 = (ny+1)*(i-1)+j;
       n2 = (ny+1)* i +j;
       edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
       u=U(edof);
   Strain_energy(ne) =u'*KE*u;
   dC(ne)=-(P*rho_old(ne)^(P-1))*Strain_energy(ne);
   V(ne)=rho_old(ne)*a*b;
   dVolume_fraction (ne)=a*b ;
   ne=ne+1;
   end
   end
   C=F'*U;
   Volume_fraction= sum(V);
   

