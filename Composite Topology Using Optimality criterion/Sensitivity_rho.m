function [C_element , dC_drho] = Sensitivity_rho(rho_1)
global Amat dof KE_A11 KE_A22 KE_A66 KE_A12 KE_A16 KE_A26 u P KE U
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[Amat,UU1,UU2,UU3,UU4,UU5] = xval_Amat(xval,Q,1,P);
A11 = Amat(1);
A12 = Amat(2);
A16 = Amat(3);
A26 = Amat(4);
A22 = Amat(5);
A66 = Amat(6);
rho = rho_1;
u = U(dof);
C_element = u'*KE*u;
dC_drho = - P * rho^-1 * u' * (A11*KE_A11+A22*KE_A22+A66*KE_A66+A12*KE_A12+A16*KE_A16+A26*KE_A26)* u;
