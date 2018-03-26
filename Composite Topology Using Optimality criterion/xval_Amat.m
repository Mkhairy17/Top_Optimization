function [Amat,U1,U2,U3,U4,U5] = xval_Amat(xval,Q,Ne,P)
U1 = 1/8 *(3*Q(1,1)+3*Q(2,2)+2*Q(1,2)+4*Q(3,3));
U2 = 1/2 *(Q(1,1)- Q(2,2));
U3 = 1/8 *(Q(1,1)+Q(2,2)-2*Q(1,2)-4*Q(3,3));
U4 = 1/8 *(Q(1,1)+Q(2,2)+6*Q(1,2)-4*Q(3,3));
U5 = 1/8 *(Q(1,1)+Q(2,2)-2*Q(1,2)+4*Q(3,3));

v1 =  xval(1:Ne,1);
v2 =  xval(Ne+1:2*Ne,1);
v3 =  xval(2*Ne+1:3*Ne,1);
v4 =  xval(3*Ne+1:4*Ne,1);
v5 =  1;
%A11,A12,A16,A26,A22,A66
Amat(:,1) = v5.^P.*(U1*ones(Ne,1) + U2.*v1 + U3.*v3);
Amat(:,2) = v5.^P.*(U4*ones(Ne,1) - U3.*v3);
Amat(:,3) = v5.^P.*(U2.*v2/2 + U3.*v4);
Amat(:,4) = v5.^P.*(U2.*v2/2 - U3.*v4);
Amat(:,5) = v5.^P.*(U1*ones(Ne,1) - U2.*v1 + U3.*v3);
Amat(:,6) = v5.^P.*(U5*ones(Ne,1)- U3.*v3);

