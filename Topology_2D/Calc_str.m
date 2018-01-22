%Calculate the strains
function strain = Calc_str(a,b,nx,ny,U)
Be = 0.5*[-1/a 0 1/a 0 1/a 0 -1/a 0 ; 0 -1/b 0 -1/b 0 1/b 0 1/b;-1/b -1/a -1/b 1/a 1/b 1/a 1/b -1/a];
strain = zeros(nx*ny,3);
ie = 1;
for ely = 1:ny
    for elx = 1:nx
      n1 = (ny+1)*(elx-1)+ely;
      n2 = (ny+1)* elx +ely;
      Ue = full(U([2*n1+1 ; 2*n1+2 ; 2*n2+1 ; 2*n2+2 ; 2*n2-1 ; 2*n2 ; 2*n1-1 ; 2*n1],1));
      strain(ie,:) = Be*Ue;
      ie = ie + 1;
    end  
end
