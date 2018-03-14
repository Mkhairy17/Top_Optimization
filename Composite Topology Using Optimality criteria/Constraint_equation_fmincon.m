 function [out1,ceq] = Constraint_equation_fmincon(xval)
out1 = zeros(2,1);      
out1(1) = 2* xval(1)^2 * (1 - xval(3)) + 2*xval(2)^2 *(1+xval(3)) + xval(3)^2 + xval(4)^2 - 4*xval(1) *xval(2) * xval(4) - 1 ;
out1(2) = xval(1)^2 + xval(2)^2 - 1 ; 
ceq = [];
      




          