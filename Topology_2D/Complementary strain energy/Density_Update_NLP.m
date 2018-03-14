function x = Density_Update_NLP(rho,nx,ny,Lx,Ly)
a = Lx/nx;
b = Ly/ny;
lb = 1e-3*ones(nx*ny,1);
ub = ones(nx*ny,1);
Aeq = [];
Beq = [];
bb = [1 ; 1]

options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
x = fmincon(@sensitivity_analysis,rho,A,bb,Aeq,Beq,lb,ub,[],options);
 