function [rho_filtered,bn,C] = Filtering_subroutine (rho_old,F,nx,ny,a,b,P)
rho_filtered = F*rho_old;
[C,dC_drho,Strain_energy] = sensitivity_analysis_optimality(rho_filtered,nx,ny,a,b,P);
dC_dx = F*dC_drho;
bn = (-dC_dx./P).*(rho_old.^(P+1));
