function [C_element,Sensitivities] = cal_C_element(xval)  
global Ne ny a b i j el P rho Amat U Ae  Lambda 
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[Amat,UU1,UU2,UU3,UU4,UU5] = xval_Amat(xval,Q,Ne);
A11 = Amat(el,1);
A12 = Amat(el,2);
A16 = Amat(el,3);
A26 = Amat(el,4);
A22 = Amat(el,5);
A66 = Amat(el,6);
Ael = [Amat(el,1) Amat(el,2) Amat(el,3);Amat(el,2) Amat(el,5) Amat(el,4);Amat(el,3) Amat(el,4) Amat(el,6) ];
[KE_A11,KE_A22,KE_A66,KE_A12,KE_A16,KE_A26] = Matrix_derivatives(a,b);
n1 = (ny+1)*(i-1)+j;
n2 = (ny+1)* i +j;
edof = [2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1]';
u = U(edof);
dC_dA11 = - rho(el)^P *u'* KE_A11 * u;
dC_dA22 = - rho(el)^P* u'* KE_A22 * u;
dC_dA66 = - rho(el)^P* u'* KE_A66 * u;
dC_dA26 = - rho(el)^P* u'* KE_A26 * u;
dC_dA16 = - rho(el)^P* u'* KE_A16 * u;
dC_dA12 = - rho(el)^P* u'* KE_A12 * u;
Psi = - [dC_dA11 0.5*dC_dA12 0.5*dC_dA16;0.5*dC_dA12 dC_dA22 0.5*dC_dA26;0.5*dC_dA16 0.5*dC_dA26 dC_dA66];
Phi = Ael * Psi* Ael;
C_element = sum(sum(Phi.*Ael^-1)) + Lambda*Ae*rho(el);


dC_drho = -P * rho(el)^(P-1) * u' * (A11*KE_A11+A22*KE_A22+A66*KE_A66+A12*KE_A12+A16*KE_A16+A26*KE_A26)* u ;
dC_dV1A = -rho(el)^(P) * u' * (UU2*KE_A11 - UU2*KE_A22)* u ;
dC_dV2A = -rho(el)^(P) * u' * (UU2/2)*(KE_A16 + KE_A26)* u ;
dC_dV3A = -rho(el)^(P) * u' * (UU3 * KE_A11 + UU3*KE_A22 - UU3*KE_A66 - UU3* KE_A12)* u ;
dC_dV4A = -rho(el)^(P) * u' * (UU3 * KE_A16 - UU3 * KE_A26)* u ;
Sensitivities = [dC_drho ; dC_dV1A ; dC_dV2A ; dC_dV3A ; dC_dV4A];      