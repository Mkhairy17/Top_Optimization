function [Strain] = Find_Strain(KE,rho_1,Lx,Ly,nx,ny,P)
a = Lx/nx;
b = Ly/ny;
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
K = global_matrix3 (KE,nx,ny,P,rho_1);
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;
Strain = Calc_str(a,b,nx,ny,U);

