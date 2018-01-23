function [rho_new,Lambda_new] = Solve_Volume_Fraction(Lambda_old,Volume_Fraction_constraint,rho_old,P,rho_min,R,a,b,nx,ny,strain,Cm,F)
rho0 = rho_old;
Volume_Fraction_old = updatedensityC(Lambda_old,rho0,R,a,b,nx,ny,P,strain,rho_min,Cm,F);
while 1
drho_dLambda = -rho_old / ((P+1)*Lambda_old);
drho_dLambda ( rho_old > 1) = 0;
drho_dLambda ( rho_old < rho_min) = 0;
deta_dLambda = sum(drho_dLambda)/(nx*ny);
Lambda_new = Lambda_old - ( Volume_Fraction_old - Volume_Fraction_constraint)/(deta_dLambda);
[Volume_Fraction_new,rho_new] = updatedensityC(Lambda_new,rho0,R,a,b,nx,ny,P,strain,rho_min,Cm,F);
if norm(Volume_Fraction_new - Volume_Fraction_old,'inf') <= 1e-10
break;
end
rho_old=rho_new;
Lambda_old=Lambda_new;
Volume_Fraction_old = Volume_Fraction_new;
end