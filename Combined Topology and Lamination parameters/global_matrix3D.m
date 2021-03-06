function [K] = global_matrix3D(Cmat,a,b,c,nx,ny,nz,P,rho)
ne = nx*ny*nz;
edofMat = DOF3D(nx,ny,nz);
row = zeros(576*ne,1);
col = zeros(576*ne,1);
val = zeros(576*ne,1);
indx = 1:576;
ie = 1;
for k = 1:nz
   for j = 1:ny
       for i = 1:nx
       [KE] = Elementstiffness_3D(Cmat(ie,:),a,b,c);
       edof = edofMat(ie,:);
       row(indx) = reshape(repmat(edof,24,1)',576,1);
       col(indx) = reshape(repmat(edof,24,1),576,1);
       val(indx) = (rho(ie)^P)*KE(:);
       indx = indx + 576;
       ie=ie+1;
       end
   end
end
 
K = sparse(row,col,val);
