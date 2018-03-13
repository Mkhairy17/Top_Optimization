function [F,rho_analysis] =Calc_F_trial(x,R,a,b,nx,ny)
% R=0.5;
% a=0.25;
% b=0.25;
% nx= 4;
% ny=4;
% x=0.5*rand(nx*ny,1);
steps_x = floor(R/a);
steps_y = floor(R/b);
element_num = 1;
x_reshaped = reshape(x,nx,ny);
x_virtual = [zeros(steps_y,nx+2*steps_x);zeros(ny,steps_x),reshape(x,nx,ny),zeros(ny,steps_x);zeros(steps_y,nx+2*steps_x)];
x_virtual=reshape(x_virtual',(ny+2*steps_y)*(nx+2*steps_x),1);
Nonzero_Indices = find(x_virtual);
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
cols = zeros(size(x_reshaped,2),1);
rows = zeros(size(x_reshaped,1),1);
indx = 1:size(Phi,1)*size(Phi,2);
for ely = 1:ny
    for elx = 1:nx
        for i = ely:2*steps_y+ely
            for j = elx:2*steps_x+elx 
                Indices = [Indices (i-1)*size(x_virtual,2)+j]; 
            end
        end
        cols(indx) = Indices;
        rows(indx) = (ely-1)*size(x_reshaped,2)+elx;
        val(indx) = Phi;
        indx = indx + size(Phi,1)*size(Phi,2); 
        Indices = [];
    end
end
F = sparse(rows,cols,val);
F = F(:,Nonzero_Indices);
F = F./sum(F,2);
rho_analysis = F*reshape(x_virtual,size(x_virtual,1)*size(x_virtual,2),1);


