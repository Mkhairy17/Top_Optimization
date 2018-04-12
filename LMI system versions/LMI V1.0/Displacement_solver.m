 function [U] = Displacement_solver(Amat,nx,ny,a,b)
 global K 
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
F(2*(nx+1)*(ny+1),1) =-1e+6;
%% Global stiffness matrix
[K]= global_matrix3(Amat,nx,ny,a,b);
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;


