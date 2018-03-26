function x = Varivbles_Update_Composite_old(xval)
lb = -1*ones(5,1);
lb(end) = 0.01;
ub = 1*ones(5,1);
Aeq = [];
Beq = [];
bb = [];
A = [];
nn = @Constraint_equation_fmincon;
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
x = fmincon(@Sensitivity_and_compliance_old,xval,A,bb,Aeq,Beq,lb,ub,nn,options);
 