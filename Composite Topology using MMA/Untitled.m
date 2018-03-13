v1 = [-1 : 0.1 : 1];
C1 = (2*v1.*v1 - sqrt(4*v1.*v1-4*(2*v1.*v1-1)))/2; %Constraint1
jbfill(v1,C1,ones(length(v1),1)')
title('Miki Diagram')


