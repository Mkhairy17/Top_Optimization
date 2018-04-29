function [C,dC] = Sensitivities_densities(U,Amat,rho,a,b,nx,ny,P)
%% Sensitivity Analysis
Strain_energy = zeros(nx*ny,1);
C = 0 ;
dC = zeros(nx*ny,1);
ne = 1;
   for j = 1:ny
   for i = 1:nx
       [KE] = Elementstiffness_composites(a,b,Amat(ne,:));
       n1 = (ny+1)*(i-1)+j;
       n2 = (ny+1)* i +j;
       edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
       u = U(edof);
       Strain_energy(ne) = u'*KE*u;
       dC(ne) = -(P*rho(ne)^(P-1))*Strain_energy(ne);
       C = C + Strain_energy(ne);
       ne = ne+1;
   end
   end
   
   
   

