function [V,Phi,L] = cal_C_element(Amat,i,j,el,ny,U,KE_A11,KE_A12,KE_A16,KE_A26,KE_A22,KE_A66,rho,P)  
%% Amat calculations
A11 = Amat(el,1);
A12 = Amat(el,2);
A16 = Amat(el,3);
A26 = Amat(el,4);
A22 = Amat(el,5);
A66 = Amat(el,6);
%%
Ael = [A11 A12 A16;A12 A22 A26 ;A16 A26 A66];
n1 = (ny+1)*(i-1)+j;
n2 = (ny+1)* i +j;
edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
u = U(edof);
dC_dA11 = -u'* KE_A11 * u;
dC_dA22 = -u'* KE_A22 * u;
dC_dA66 = -u'* KE_A66 * u;
dC_dA26 = -u'* KE_A26 * u;
dC_dA16 = -u'* KE_A16 * u;
dC_dA12 = -u'* KE_A12 * u;
Psi = - [dC_dA11 0.5*dC_dA12 0.5*dC_dA16;0.5*dC_dA12 dC_dA22 0.5*dC_dA26;0.5*dC_dA16 0.5*dC_dA26 dC_dA66];
Phi = Ael * Psi* Ael * rho^P;
[V,D] = eig(Phi);
V = V*sqrt(D);
L = chol(Phi);