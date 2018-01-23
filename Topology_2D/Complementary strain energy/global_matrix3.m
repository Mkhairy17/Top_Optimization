function [K]=global_matrix3(KE,nx,ny,P,rho)
ne = nx*ny;
row = zeros(64*ne,1);
col = zeros(64*ne,1);
val = zeros(64*ne,1);
indx = 1:64;
ie=1;
for j = 1:ny
   for i = 1:nx
       n1 = (ny+1)*(i-1)+j;
       n2 = (ny+1)* i +j;
       edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
       row(indx) = reshape(repmat(edof,8,1)',64,1);
       col(indx) = reshape(repmat(edof,8,1),64,1);
       val(indx) = (rho(ie)^P)*KE(:);
       indx = indx + 64;
       ie=ie+1;
    end
end    
 
K = sparse(row,col,val);
