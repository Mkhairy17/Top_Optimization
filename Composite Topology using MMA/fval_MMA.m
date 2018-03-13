 function [out1] = fval_MMA(xval,eita_V0,Ae,Ne)

out1 = zeros(2*Ne+1 ,1);      
for i = 1:Ne
    out1(2*i-1) = 2* xval(i)^2 * (1 - xval(i+2*Ne)) + 2*xval(i+Ne)^2 *(1+xval(i+2*Ne)) + xval(i+2*Ne)^2 + xval(i+3*Ne)^2 - 4*xval(i) *xval(i+Ne) * xval(i+3*Ne) - 1 ;
    out1(2*i) = xval(i)^2 + xval(i+Ne)^2 - 1 ; 
end
out1(2*Ne+1) = sum(xval(4*Ne+1:5*Ne)) * Ae - eita_V0;
      




          