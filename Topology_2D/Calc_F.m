function [F] = Calc_F(rho_1,strain,R,a,b,nx,ny,P,Cm)
ne = size(strain,1);
steps_x = floor(R/a);
steps_y = floor(R/b);
for ie = 1:ne
   Strain_energy(ie) = rho_1(ie)^(2*P)*(strain(ie,:)*Cm*strain(ie,:)'*a*b); 
end
Strain_energy_reshaped = reshape(Strain_energy,nx,ny)';
Strain_energy_reshaped_template = [zeros(steps_y,nx+2*steps_x);zeros(ny,steps_x),Strain_energy_reshaped,zeros(ny,steps_x);zeros(steps_y,nx+2*steps_x)];
element_num = 1;
for r_x = -a*steps_x:a:a*steps_x
    for r_y = -b*steps_y:b:b*steps_y  
        r = sqrt(r_x^2 + r_y^2);
        Phi(element_num) = 1 - (r/R)^2;
        element_num = element_num + 1;
    end
end
Phi(Phi < 0) = 0;
Phi = Phi/sum(sum(Phi));
Indices = [];
cols = zeros(size(Strain_energy_reshaped,2),1);
rows = zeros(size(Strain_energy_reshaped,1),1);
indx = 1:size(Phi,1)*size(Phi,2);
for ely = 1:ny
    for elx = 1:nx
        for i = ely:2*steps_y+ely
            for j = elx:2*steps_x+elx 
                Indices = [Indices (i-1)*size(Strain_energy_reshaped_template,2)+j]; 
            end
        end
        cols(indx) = Indices;
        rows(indx) = (ely-1)*size(Strain_energy_reshaped,2)+elx;
        val(indx) = Phi;
        indx = indx + size(Phi,1)*size(Phi,2); 
        Indices = [];
    end
end
F = sparse(rows,cols,val);
