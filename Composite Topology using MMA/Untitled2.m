Q = [181.81 2.9 0; 2.9 10.35  0;0 0 7.17] * 1e9;
U1 = 1/8 *(3*Q(1,1)+3*Q(2,2)+2*Q(1,2)+4*Q(3,3));
U2 = 1/2 *(Q(1,1)- Q(2,2));
U3 = 1/8 *(Q(1,1)+Q(2,2)-2*Q(1,2)-4*Q(3,3));
U4 = 1/8 *(Q(1,1)+Q(2,2)+6*Q(1,2)-4*Q(3,3));
U5 = 1/8 *(Q(1,1)+Q(2,2)-2*Q(1,2)+4*Q(3,3));
C2 = ((U4^2-U1^2) + U2.^2.*v1.^2)./((2*U3)*(U1+U4)).*ones(length(v1),1)';
C3 = -(U1+U2*v1)./(U3*ones(length(v1),1)');
C4 = U5/U3*ones(length(v1),1)';
jbfill(v1,C2,C4);
hold on
jbfill(v1,C3,C4);
title('Positive definivie matrix constaints')