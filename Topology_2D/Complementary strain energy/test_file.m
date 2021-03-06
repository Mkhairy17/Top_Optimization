clear; clc
global nx ny Lx Ly P 
nx=40;    ny=40;
Lx=1;     Ly=1;
P=1;       rho = 0.5*ones(nx,ny);
[C,dC,Volume_fraction, dVolume_fraction] = sensitivity_analysis(rho,Lx,Ly,nx,ny,P);
A=dVolume_fraction;   b= 0.5;    lb= 1e-3*ones(nx*ny,1);   ub=1*ones(nx*ny,1);    Aeq=[];     Beq=[];         x0=rho;  nonlcon = [];
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
fun = @sensitivity_analysis2;
x = fmincon(fun,x0,A,b,Aeq,Beq,lb,ub,nonlcon,options);
colormap(gray); imagesc(-x'); axis equal; axis tight; axis off;pause(1e-6);
 