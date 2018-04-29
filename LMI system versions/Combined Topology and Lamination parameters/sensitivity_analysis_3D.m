function [dC] = sensitivity_analysis_3D(U,Cmat,rho,nx,ny,nz,a,b,c,P)
%% Sensitivity Analysis
edofMat = DOF3D(nx,ny,nz);
Strain_energy = zeros(nx*ny*nz,1);
Strain_energy_penalized = zeros(nx*ny*nz,1);
dC = zeros(nx*ny*nz,1);
ne = 1;
for j = 1:nz
    for i = 1:ny
        for k = 1:nx
       [KE] = Elementstiffness_3D(Cmat(ne,:),a,b,c);
       edof = edofMat(ne,:);
       u = U(edof);
       Strain_energy(ne) = u'*KE*u;
       Strain_energy_penalized(ne) = rho(ne)^P*u'*KE*u;
       dC(ne) = -(P*rho(ne)^(P-1))*Strain_energy(ne);
       ne = ne+1;
        end
    end
end
  
   
   

