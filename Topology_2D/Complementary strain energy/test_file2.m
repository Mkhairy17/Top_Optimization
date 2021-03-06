function x = Density_Update_NLP(rho,nx,ny,Lx,Ly,P)


a = Lx/nx;
b = Ly/ny;
lb = 1e-3*ones(nx*ny,1);
ub = ones(nx*ny,1);
Aeq = [];
Beq = [];

bb = 0.5;
A = a*b*ones(1,nx*ny);
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
x = fmincon(@sensitivity_analysis,rho,A,bb,Aeq,Beq,lb,ub,[],options);
 