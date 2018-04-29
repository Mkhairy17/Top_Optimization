function [rho_filtered,bn] = Filtering_subroutine_3D(U,Cmat,rho_old,F,nx,ny,nz,a,b,c,P)
rho_filtered = F*rho_old;
[dC_drho] = sensitivity_analysis_3D(U,Cmat,rho_old,nx,ny,nz,a,b,c,P);
dC_dx = F*dC_drho;
bn = (-dC_dx./P).*(rho_old.^(P+1));
