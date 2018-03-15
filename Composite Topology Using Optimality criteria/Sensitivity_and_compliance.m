function [C_element , Sensitivities_current] = Sensitivity_and_compliance(xval)
global Sensitivities Amat Phi rho Lambda Ne el Ae
Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
[Amat,UU1,UU2,UU3,UU4,UU5] = xval_Amat(xval,Q,1);
A11 = Amat(1);
A12 = Amat(2);
A16 = Amat(3);
A26 = Amat(4);
A22 = Amat(5);
A66 = Amat(6);
Ael = [A11 A12 A16;A12 A22 A26 ;A16 A26 A66 ];C_element = sum(sum(Phi.*Ael^-1)) + Lambda*Ae*rho(el);
Sensitivities_current = Sensitivities;