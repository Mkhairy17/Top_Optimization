function [rho_filtered,bn,C] = Filtering_subroutine (Amat,U,rho_old,F,nx,ny,a,b,P)
rho_filtered = F*rho_old;
[C,dC_drho] = Sensitivities_densities(U,Amat,rho_old,a,b,nx,ny,P);
dC_dx = F*dC_drho;
bn = (-dC_dx./P).*(rho_old.^(P+1));
