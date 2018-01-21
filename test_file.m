clear; clc
nx=5;    ny=5;
Lx=1;     Ly=1;
P=1;       rho = 0.5*ones(nx,ny);
[C,dC,Volume_fraction, dVolume_fraction] = sensitivity_analysis(rho,Lx,Ly,nx,ny,P);
 