clear
clc
nx = 4; 
ny = 4; 
nz = 4;
ne = 1; 
Ne = nx*ny*nz;
parfor j = 1:ny
    for i = 1:nx
        for k = 1:nz
%         edof = edofMat(ne,:);
%         [L,phi] = cal_C_element(C,edof,U,KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66,rho_old(ne),P);
        a1 = 1;
        a2 = 1; 
        a3 = 1;
        a4 = 1;
        a5 = 1;
        L = i*j*k* ones(6,6);
        xval_opt = optimum_variables (a1,a2,a3,a4,a5,L);
%         xval_new(i) = xval_opt(4);
%         xval_new(ne+Ne) = xval_opt(5);
%         xval_new(ne+2*Ne) = xval_opt(6);
%         xval_new(ne+3*Ne) = xval_opt(7); 
%         i 
%         j 
%         k
%         ne
%         i+j+k-2
%         ne = ne+1;
        end
    end
end