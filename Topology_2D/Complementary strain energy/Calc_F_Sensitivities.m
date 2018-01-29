function [F] = Calc_F_Sensitivities(Filter_template,R,a,b,nx,ny,P)
steps_x = floor(R/a);
steps_y = floor(R/b);
Filter_template_reshaped = reshape(Filter_template,nx,ny)';
Strain_energy_reshaped_template = [zeros(steps_y,nx+2*steps_x);zeros(ny,steps_x),Filter_template_reshaped,zeros(ny,steps_x);zeros(steps_y,nx+2*steps_x)];
Filter_template_reshaped_template_vec = reshape(Strain_energy_reshaped_template',(ny+2*steps_y)*(nx+2*steps_x),1);
Nonzero_Indices = find(Filter_template_reshaped_template_vec);
element_num = 1;
for r_x = -a*steps_x:a:a*steps_x
    for r_y = -b*steps_y:b:b*steps_y  
        r = sqrt(r_x^2 + r_y^2);
        Phi(element_num) = 1 - (r/R)^2;
        element_num = element_num + 1;
    end
end
Phi(Phi < 0) = 0;
% Phi = Phi/sum(sum(Phi));
Indices = [];
cols = zeros(size(Filter_template_reshaped,2),1);
rows = zeros(size(Filter_template_reshaped,1),1);
indx = 1:size(Phi,1)*size(Phi,2);
for ely = 1:ny
    for elx = 1:nx
        for i = ely:2*steps_y+ely
            for j = elx:2*steps_x+elx 
                Indices = [Indices (i-1)*size(Strain_energy_reshaped_template,2)+j]; 
            end
        end
        cols(indx) = Indices;
        rows(indx) = (ely-1)*size(Filter_template_reshaped,2)+elx;
        val(indx) = Phi;
        indx = indx + size(Phi,1)*size(Phi,2); 
        Indices = [];
    end
end
F1 = sparse(rows,cols,val);
F = F1(:,Nonzero_Indices);
F = F./sum(F,2);
