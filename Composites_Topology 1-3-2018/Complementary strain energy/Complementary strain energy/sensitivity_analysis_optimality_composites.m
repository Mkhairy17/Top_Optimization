function [C,dC_drho,dC_dV1A,dC_dV2A,dC_dV3A,dC_dV4A,Strain_energy] = sensitivity_analysis_optimality_composites(rho,Amat,UU2,UU3,nx,ny,a,b,P)
%% Parameters
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
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
[K]= global_matrix3(Amat,rho,nx,ny,P,a,b);
%Define displacement vector
U = sparse(2*(ny+1)*(nx+1),1);
%Solve for U
U(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
U(FixDOF,:) = 0;

%% Sensitivity Analysis
Strain_energy = zeros(nx*ny,1);
dC_dV1A = zeros(nx*ny,1);
dC_dV2A = zeros(nx*ny,1);
dC_dV3A = zeros(nx*ny,1);
dC_dV4A = zeros(nx*ny,1);
dC_drho = zeros(nx*ny,1);
ne = 1;
   for j = 1:ny
   for i = 1:nx
       KE = Elementstiffness_composites(a,b,Amat(ne,:)); %A11,A12,A16,A26,A22,A66
       A11 = Amat(ne,1);
       A12 = Amat(ne,2);
       A16 = Amat(ne,3);
       A26 = Amat(ne,4);
       A22 = Amat(ne,5);
       A66 = Amat(ne,6);
       n1 = (ny+1)*(i-1)+j;
       n2 = (ny+1)* i +j;
       edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
       u = U(edof);
       Strain_energy(ne) = u'*KE*u;
       dC_drho(ne) = -P * rho(ne)^(2*P-1) * u' * (A11*KE_A11+A22*KE_A22+A66*KE_A66+A12*KE_A12+A16*KE_A16+A26*KE_A26)* u ;
       dC_dV1A(ne) = -rho(ne)^(P) * u' * (UU2*KE_A11 - UU2*KE_A22)* u ;
       dC_dV2A(ne) = -rho(ne)^(P) * u' * (UU2/2)*(KE_A16 + KE_A26)* u ;
       dC_dV3A(ne) = -rho(ne)^(P) * u' * (UU3 * KE_A11 + UU3*KE_A22 - UU3*KE_A66 - UU3* KE_A12)* u ;
       dC_dV4A(ne) = -rho(ne)^(P) * u' * (UU3 * KE_A16 - UU3 * KE_A26)* u ;
       ne = ne+1;
   end
   end
   C = F'*U;