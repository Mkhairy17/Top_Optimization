function [Strain_energy_new,Avg_Lambda0] = Calc_Strain_Energy(rho_1,strain,R,P,nx,ny,a,b,Cm,F)
ne = size(strain,1);
steps_x = floor(R/a);
steps_y = floor(R/b);
for ie = 1:ne
   Strain_energy(ie) = rho_1(ie)^(2*P)*(strain(ie,:)*Cm*strain(ie,:)'*a*b); 
end
Strain_energy_reshaped = reshape(Strain_energy,nx,ny)';
Strain_energy_reshaped_template = [zeros(steps_y,nx+2*steps_x);zeros(ny,steps_x),Strain_energy_reshaped,zeros(ny,steps_x);zeros(steps_y,nx+2*steps_x)];
Strain_energy_reshaped_template_col = reshape(Strain_energy_reshaped_template',1,(ny+2*steps_y)*(nx+2*steps_x));
Strain_energy_new = F*Strain_energy_reshaped_template_col';
