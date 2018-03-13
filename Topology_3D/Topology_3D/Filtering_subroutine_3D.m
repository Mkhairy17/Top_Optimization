function [rho_filtered,bn,C] = Filtering_subroutine_3D(rho_old,F,nx,ny,nz,Lx,Ly,Lz,P)
rho_filtered = F*rho_old;
[C,dC_drho] = sensitivity_analysis_3D(rho_filtered,nx,ny,nz,Lx,Ly,Lz,P);
dC_dx = F*dC_drho;
bn = (-dC_dx./P).*(rho_old.^(P+1));
