function [L,Phi] = cal_C_element(C,edof,U,KE_C11,KE_C12,KE_C13,KE_C14,KE_C15,KE_C16,KE_C22,KE_C23,KE_C24,KE_C25,KE_C26,KE_C33,KE_C34,KE_C35,KE_C36,KE_C44,KE_C45,KE_C46,KE_C55,KE_C56,KE_C66,rho_element,el,P) 
%% C calculations
C11 = C(el,1);
C12 = C(el,2);
C13 = C(el,3);
C14 = C(el,4);
C15 = C(el,5);
C16 = C(el,6);
C22 = C(el,7);
C23 = C(el,8);
C24 = C(el,9);
C25 = C(el,10);
C26 = C(el,11);
C33 = C(el,12);
C34 = C(el,13);
C35 = C(el,14);
C36 = C(el,15);
C44 = C(el,16);
C45 = C(el,17);
C46 = C(el,18);
C55 = C(el,19);
C56 = C(el,20);
C66 = C(el,21);

%%
CC = [C11 C12 C13 C14 C15 C16 ;C13 C22 C23 C24 C25 C26;C13 C23 C33 C34 C35 C36;C14 C24 C34 C44 C45 C46;C15 C25 C35 C45 C55 C56;C16 C26 C36 C46 C56 C66];
u = U(edof);
dC_dC11 = -u'* KE_C11 * u;
dC_dC12 = -u'* KE_C12 * u;
dC_dC13 = -u'* KE_C13 * u;
dC_dC15 = -u'* KE_C15 * u;
dC_dC16 = -u'* KE_C16 * u;
dC_dC22 = -u'* KE_C22 * u;
dC_dC23 = -u'* KE_C23 * u;
dC_dC25 = -u'* KE_C25 * u;
dC_dC26 = -u'* KE_C26 * u;
dC_dC33 = -u'* KE_C33 * u;
dC_dC34 = -u'* KE_C34 * u;
dC_dC35 = -u'* KE_C35 * u;
dC_dC36 = -u'* KE_C36 * u;
dC_dC14 = -u'* KE_C14 * u;
dC_dC24 = -u'* KE_C24 * u;
dC_dC44 = -u'* KE_C44 * u;
dC_dC45 = -u'* KE_C45 * u;
dC_dC46 = -u'* KE_C46 * u;
dC_dC55 = -u'* KE_C55 * u;
dC_dC56 = -u'* KE_C56 * u;
dC_dC66 = -u'* KE_C66 * u;
%%=======**=========


Psi = - [dC_dC11 0.5*dC_dC12 0.5*dC_dC13 0.5*dC_dC14 0.5*dC_dC15 0.5*dC_dC16; 0.5*dC_dC12 dC_dC22 0.5*dC_dC23 0.5*dC_dC24 0.5*dC_dC25 0.5*dC_dC26 ;0.5*dC_dC13 0.5*dC_dC23 dC_dC33 0.5*dC_dC34 0.5*dC_dC35 0.5*dC_dC36;0.5*dC_dC14 0.5*dC_dC24 0.5*dC_dC34 dC_dC44 0.5*dC_dC45 0.5*dC_dC46;0.5*dC_dC15 0.5*dC_dC25 0.5*dC_dC35 0.5*dC_dC45 dC_dC55 0.5*dC_dC56;0.5*dC_dC16 0.5*dC_dC26 0.5*dC_dC36 0.5*dC_dC46 0.5*dC_dC56 dC_dC66];
Phi = CC * Psi* CC *rho_element^P;
L = real(sqrtm(Phi));
% [V,D] = eig(Phi);
% V = V*sqrt(D);