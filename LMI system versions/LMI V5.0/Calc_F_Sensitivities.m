function [F] = Calc_F_Sensitivities(Filter_template,R,a,b,nx,ny)
steps_x = floor(R/a);
steps_y = floor(R/b);
Filter_template_reshaped = reshape(Filter_template,nx,ny)';
New_template_reshaped = [zeros(steps_y,nx+2*steps_x);zeros(ny,steps_x),Filter_template_reshaped,zeros(ny,steps_x);zeros(steps_y,nx+2*steps_x)];
Filter_template_reshaped_vec = reshape(New_template_reshaped',(ny+2*steps_y)*(nx+2*steps_x),1);
Nonzero_Indices = find(Filter_template_reshaped_vec);
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
cols = zeros(length(Phi)*nx*ny,1);
rows = zeros(length(Phi)*nx*ny,1);
indx = 1:size(Phi,1)*size(Phi,2);
for ely = 1:ny
    for elx = 1:nx
        for j = ely:2*steps_y+ely
            for i = elx:2*steps_x+elx 
                Indices = [Indices (j-1)*(nx+2*steps_x)+i]; 
            end
        end
        cols(indx) = Indices;
        rows(indx) = (ely-1)*nx+elx;
        val(indx) = Phi;
        indx = indx + size(Phi,1)*size(Phi,2); 
        Indices = [];
    end
end
F1 = sparse(rows,cols,val);
F = F1(:,Nonzero_Indices);
F = F./sum(F,2);
