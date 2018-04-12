function [C,dC] = sensitivity_analysis_3D(rho,nx,ny,nz,Lx,Ly,Lz,P)
a = Lx/nx; %elemeny width
b = Ly/ny; %element height
c = Lz/nz;  %Element depth
%% Parameters
E = 1;
v = 0.3;
KE = Elementstiffness_3D2(a,b,c,E,v);
K = global_matrix3D(KE,nx,ny,nz,P,rho);
ndof = 3*(nx+1)*(ny+1)*(nz+1);
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nx, 0, 0:nz);                 % Coordinates
loadnid = kl*(nx+1)*(ny+1)+il*(ny+1)+(ny+1-jl);     % Node IDs
loaddof = 3*loadnid(:) - 1;                         % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:ny,0:nz);                  % Coordinates
fixednid = kf*(nx+1)*(ny+1)+iif*(ny+1)+(ny+1-jf);     % Node IDs
FixDOF = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
%% Force vector
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
FreeDOF = setdiff(1:ndof,FixDOF);
%Solve for U
% U_3D(FreeDOF,:) = K(FreeDOF,FreeDOF) \ F(FreeDOF,:);
% U_3D(FixDOF,:) = 0;
tolit = 1e-8;
maxit = 8000;
M = diag(diag(K(FreeDOF,FreeDOF))) ;
U(FreeDOF,:) = pcg (K(FreeDOF,FreeDOF),F(FreeDOF,:),tolit,maxit,M) ;


%% Sensitivity Analysis
edofMat = DOF3D(nx,ny,nz);
Strain_energy = zeros(nx*ny*nz,1);
Strain_energy_penalized = zeros(nx*ny*nz,1);
dC = zeros(nx*ny*nz,1);
ne = 1;
for j = 1:nz
    for i = 1:ny
        for k = 1:nx
       edof = edofMat(ne,:);
       u = U(edof);
       Strain_energy(ne) = u'*KE*u;
       Strain_energy_penalized(ne) = rho(ne)^P*u'*KE*u;
       dC(ne) = -(P*rho(ne)^(P-1))*Strain_energy(ne);
       ne = ne+1;
        end
    end
end
   C = F'*U;
   
   

