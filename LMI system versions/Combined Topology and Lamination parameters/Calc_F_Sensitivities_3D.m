function [F] = Calc_F_Sensitivities_3D(Filter_template,R,a,b,c,nx,ny,nz)
steps_x = floor(R/a);
steps_y = floor(R/b);
steps_z = floor(R/c);
Filter_initial = reshape(Filter_template,nx,ny,nz);
Filter_reshaped_template = zeros(nx+2*steps_x,ny+2*steps_y,nz+2*steps_z);
Filter_reshaped_template(1+steps_x : nx+steps_x,1+steps_y : ny+steps_y,1+steps_z : nz+steps_z) =  Filter_initial;
Filter_template_vec = reshape(Filter_reshaped_template,(ny+2*steps_y)*(nx+2*steps_x)*(nz+2*steps_z),1);
Nonzero_Indices = find(Filter_template_vec);
element_num = 1;
for r_x = -a*steps_x:a:a*steps_x
    for r_y = -b*steps_y:b:b*steps_y  
        for r_z = -c*steps_z:c:c*steps_z
        r = sqrt(r_x^2 + r_y^2 + r_z^2);
        Phi(element_num) = 1 - (r/R)^2;
        element_num = element_num + 1;
        end
    end
end
Phi(Phi < 0) = 0;
Indices = [];
cols = zeros(length(Phi)*nx*ny*nz,1);
rows = zeros(length(Phi)*nx*ny*nz,1);
indx = 1:size(Phi,1)*size(Phi,2);
for elz = 1:nz
    for ely = 1:ny
        for elx = 1:nx
            for k = elz:2*steps_z+elz
                for  j = ely:2*steps_y+ely
                    for  i = elx:2*steps_x+elx 
                Indices = [Indices (k-1)*(nx+2*steps_x)*(ny+2*steps_y)+(j-1)*(nx+2*steps_x)+i];
                    end
                end
            end
        cols(indx) = Indices;
        rows(indx) = (elz-1)*nx*ny+(ely-1)*nx+elx;
        val(indx) = Phi;
        indx = indx + size(Phi,1)*size(Phi,2); 
        Indices = [];
        end
    end
end
F1 = sparse(rows,cols,val);
F = F1(:,Nonzero_Indices);
F = F./sum(F,2);
