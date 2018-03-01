function [dfdx]= dfdx_calc(xval,Ne,Ae)

v1 =  xval(1:Ne,1);
v2 =  xval(Ne+1:2*Ne,1);
v3 =  xval(2*Ne+1:3*Ne,1);
v4 = xval(3*Ne+1:4*Ne,1);
df1dvi=[4*v1'-4*v1'.*v3'-4*v2'.*v4' 4*v2'+4*v2'.*v3'-4*v1'.*v4' -2*v1'.*v1'+2*v2'.*v2'+2*v3' 2*v4'-4*v1'.*v2' zeros(1,Ne)];
df2dvi=[2*v1' 2*v2' zeros(1,Ne) zeros(1,Ne) zeros(1,Ne)];     
df3dvi=[zeros(1,Ne) zeros(1,Ne) zeros(1,Ne) zeros(1,Ne) Ae'];

m=2*Ne+1;
n=5*Ne;

dfdx=zeros(m,n);

dfdx(1:2:2*Ne,:)= repmat(df1dvi,[2,1]);
dfdx(2:2:2*Ne,:)= repmat(df2dvi,[2,1]);
dfdx(end,:)=df3dvi;

    
    
    
    
        