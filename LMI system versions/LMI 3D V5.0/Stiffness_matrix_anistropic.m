%fungamation [ke] = Elementstiffness(Lex,Ley,Lez)
clc
clear 
syms eta zeta gama a b c C1111 C1122 C1133 C1123 C1113 C1112 C2222 C2233 C2223 C2213 C2212 C3333 C3323 C3313 C3312 C2323 C2313 C2312 C1313 C1312 C1212


nq_ezetasilon=0.125*[(1-eta*2/a)*(1-zeta*2/b)*(1-gama*2/c);(1+eta*2/a)*(1-zeta*2/b)*(1-gama*2/c);(1+eta*2/a)*(1+zeta*2/b)*(1-gama*2/c);(1-eta*2/a)*(1+zeta*2/b)*(1-gama*2/c);(1-eta*2/a)*(1-zeta*2/b)*(1+gama*2/c);(1+eta*2/a)*(1-zeta*2/b)*(1+gama*2/c);(1+eta*2/a)*(1+zeta*2/b)*(1+gama*2/c);(1-eta*2/a)*(1+zeta*2/b)*(1+gama*2/c)];
zeta11=diff(nq_ezetasilon(1),eta);
zeta12=diff(nq_ezetasilon(1),zeta);
zeta13=diff(nq_ezetasilon(1),gama);
zeta21=diff(nq_ezetasilon(2),eta);
zeta22=diff(nq_ezetasilon(2),zeta);
zeta23=diff(nq_ezetasilon(2),gama);
zeta31=diff(nq_ezetasilon(3),eta);
zeta32=diff(nq_ezetasilon(3),zeta);
zeta33=diff(nq_ezetasilon(3),gama);
zeta41=diff(nq_ezetasilon(4),eta);
zeta42=diff(nq_ezetasilon(4),zeta);
zeta43=diff(nq_ezetasilon(4),gama);
zeta51=diff(nq_ezetasilon(5),eta);
zeta52=diff(nq_ezetasilon(5),zeta);
zeta53=diff(nq_ezetasilon(5),gama);
zeta61=diff(nq_ezetasilon(6),eta);
zeta62=diff(nq_ezetasilon(6),zeta);
zeta63=diff(nq_ezetasilon(6),gama);
zeta71=diff(nq_ezetasilon(7),eta);
zeta72=diff(nq_ezetasilon(7),zeta);
zeta73=diff(nq_ezetasilon(7),gama);
zeta81=diff(nq_ezetasilon(8),eta);
zeta82=diff(nq_ezetasilon(8),zeta);
zeta83=diff(nq_ezetasilon(8),gama);

B=[zeta11 0 0 zeta21 0 0 zeta31 0 0 zeta41 0 0 zeta51 0 0 zeta61 0 0 zeta71 0 0 zeta81 0 0;...
    0 zeta12 0 0 zeta22 0 0 zeta32 0 0 zeta42 0 0 zeta52 0 0 zeta62 0 0 zeta72 0 0 zeta82 0;...
    0 0 zeta13 0 0 zeta23 0 0 zeta33 0 0 zeta43 0 0 zeta53 0 0 zeta63 0 0 zeta73 0 0 zeta83;...
    zeta12 zeta11 0 zeta22 zeta21 0 zeta32 zeta31 0 zeta42 zeta41 0 zeta52 zeta51 0 zeta62 zeta61 0 zeta72 zeta71 0 zeta82 zeta81 0;...
    0 zeta13 zeta12 0 zeta23 zeta22 0 zeta33 zeta32 0 zeta43 zeta42 0 zeta53 zeta52 0 zeta63 zeta62 0 zeta73 zeta72 0 zeta83 zeta82;...
    zeta13 0 zeta11 zeta23 0 zeta21 zeta33 0 zeta31 zeta43 0 zeta41 zeta53 0 zeta51 zeta63 0 zeta61 zeta73 0 zeta71 zeta83 0 zeta81];

C = [C1111 C1122 C1133 C1123 C1113 C1112 ;C1122 C2222 C2233 C2223 C2213 C2212;C1133 C2233 C3333 C3323 C3313 C3312;C1123 C2223 C3323 C2323 C2313 C2312;C1113 C2213 C3313 C2313 C1313 C1312;C1112 C2212 C3312 C2312 C1312 C1212];

ke=int(B'*C*B,eta,-a/2,a/2);
ke=int(ke,zeta,-b/2,b/2);
ke=int(ke,gama,-c/2,c/2);

 