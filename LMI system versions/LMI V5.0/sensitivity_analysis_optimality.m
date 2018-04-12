function [C,dC,Strain_energy] = sensitivity_analysis_optimality(U,rho,Amat,nx,ny,a,b,P)
%% Sensitivity Analysis
Strain_energy = zeros(nx*ny,1);
dC = zeros(nx*ny,1);
ne = 1;
   for j = 1:ny
   for i = 1:nx
       KE = Elementstiffness_composites(a,b,Amat(ne,:)); %A11,A12,A16,A26,A22,A66
       n1 = (ny+1)*(i-1)+j;
       n2 = (ny+1)* i +j;
       edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
       u = U(edof);
       Strain_energy(ne) = u'*KE*u;
       dC(ne) = -(P*rho(ne)^(P-1))*Strain_energy(ne);
       ne = ne+1;
   end
   end
   
   

