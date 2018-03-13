function [K]=global_matrix2(Lx,Ly,nx,ny,P,rho)

% dx=Lx/nx;                    %element size in x
% dy=Ly/ny;                    %element size in y 
X=nx+1;                       %number of nodes in x
Y=ny+1;                       %number of nodes in y
nnode = (nx+1)*(ny+1);        %total number of nodes in system
K = sparse(2*X*Y,2*X*Y);      %sparse matrix of the global stiffness
%local stiffness matrix
E = 1;
nu = 0.3;
k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8); k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3); k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2); k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5);k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4);k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7);k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6);k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
for j = 1:ny
   for i = 1:nx
 n1 = (ny+1)*(i-1)+j;
 n2 = (ny+1)* i +j;
 edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1;2*n2+2;2*n1+1; 2*n1+2];
 K(edof,edof) = K(edof,edof)+KE;
    end
end
