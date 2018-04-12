function U = Displacement_Solver(Cmat,rho,a,b,c,nx,ny,nz,P)
[K] = global_matrix3D(Cmat,a,b,c,nx,ny,nz,P,rho);
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