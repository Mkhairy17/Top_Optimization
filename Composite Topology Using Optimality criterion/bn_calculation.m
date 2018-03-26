function [bn] = bn_calculation (U,Amat,rho,a,b,nx,ny,P)
[C,dC] = Sensitivities_densities(U,Amat,rho,a,b,nx,ny,P);
bn = (-dC./P).*(rho.^(P+1));
