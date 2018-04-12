 function [U] = Displacement_solver_density(Amat,rho_old,nx,ny,a,b,P)
%% Parameters
% All degrees of freedom
AllDOF = 1:2*(nx+1)*(ny+1);
%Set Fixed degrees of freedom
FixDOF = [1:2*(ny+1)];
%Set free degrees of freedom
FreeDOF = setdiff(AllDOF,FixDOF);
%Define force vector
F = sparse(2*(ny+1)*(nx+1),1); 
%%Define forces
F(2*(nx+1)*(ny+1),1) = -1*10^6;
%% Global stiffness matrix
[K] = global_matrix_penalized(Amat,rho_old,nx,ny,a,b,P); %nonpenalized stiffness matrix
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;


