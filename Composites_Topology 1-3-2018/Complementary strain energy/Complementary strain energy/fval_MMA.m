function [out1]=fval_MMA(xval,eita_V0,Ae,Ne)
Matrix = [2*ones(Ne,Ne) -2*ones(Ne,Ne) 2*ones(Ne,Ne) 2*ones(Ne,Ne) Ne*ones(Ne,Ne) Ne*ones(Ne,Ne) -4*ones(Ne,Ne) zeros(Ne,Ne);
          Ne*ones(Ne,Ne)  zeros(Ne,Ne)  Ne*ones(Ne,Ne) zeros(Ne,Ne)  zeros(Ne,Ne)  zeros(Ne,Ne)   zeros(Ne,Ne)  zeros(Ne,Ne);
          zeros(1,Ne)   zeros(1,Ne)  zeros(1,Ne)  zeros(1,Ne)  zeros(1,Ne)  zeros(1,Ne)   zeros(1,Ne)  Ae*ones(1,Ne)];
      
Vector = [xval(1:Ne,1).*xval(1:Ne,1);
          xval(1:Ne,1).*xval(1:Ne,1).*xval(2*Ne+1:3*Ne,1);
          xval(Ne+1:2*Ne,1).*xval(Ne+1:2*Ne,1);
          xval(Ne+1:2*Ne,1).*xval(Ne+1:2*Ne,1).*v3;
          xval(2*Ne+1:3*Ne,1).*xval(2*Ne+1:3*Ne,1);
          xval(3*Ne+1:4*Ne,1).*xval(3*Ne+1:4*Ne,1);
          xval(1:Ne,1).*xval(Ne+1:2*Ne,1).*xval(3*Ne+1:4*Ne,1);
          xval(4*Ne+1:5*Ne,1)]; 
      
Cvector = [ones(Ne,1);ones(Ne,1);eita_V0]; 


out1 = Matrix*Vector-Cvector;


          