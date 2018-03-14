function x = Varivbles_Update_Composite(xval)
lb = -1*ones(5,1);
lb(end) = 1e-3;
ub = 1*ones(5,1);
Aeq = [];
Beq = [];
bb = [];
A = [];
nn = @Constraint_equation_fmincon;
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
x = fmincon(@cal_C_element,xval,A,bb,Aeq,Beq,lb,ub,nn,options);
 