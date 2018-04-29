function [C] = Elasticity_Matrix(a1,a2,a3,a4,a5,F1111,F1122,F1133,F2222,F2233,F3333,F2333,F2223,F1123,F1113,F1233,F1333,F1112,F1222,F1223)
C11 = 2*a1+2*a2+2*F1111*a3+2*(F1111+F1122+F1133)*a4+2*(F1111+F1122+F1133)*a5; 
C12 = 2*a1+2*F1122*a3+(F1111+F1122+F1133)*a4+(F1122+F2222+F2233)*a4;
C13 = 2*a1+2*F1133*a3+(F1111+F1122+F1133)*a4+(F1133+F2233+F3333)*a4;
C14 = 4*F1123*a3+2*(F1123+F2223+F2333)*a4;
C15 = 4*F1113*a3+2*(F1113+F1233+F1333)*a4;
C16 = 4*F1112*a3+2*(F1112+F1222+F1233)*a4;
C22 = 2*a1+2*a2+2*F2222*a3+2*(F1122+F2222+F2233)*a4+2*(F1122+F2222+F2233)*a5;
C23 = 2*a1+2*F2233*a3+(F1122+F2222+F2233)*a4+(F1133+F2233+F3333)*a4;
C24 = 4*F2223*a3+2*(F1123+F2223+F2333)*a4;
C25 = 4*F1223*a3+2*(F1113+F1233+F1333)*a4;
C26 = 4*F1222*a3+2*(F1112+F1222+F1233)*a4;
C33 = 2*a1+2*a2+2*F3333*a3+2*(F1133+F2233+F3333)*a4+2*(F1133+F2233+F3333)*a5;
C34 = 4*F2333*a3+2*(F1123+F2223+F2333)*a4;
C35 = 4*F1333*a3+2*(F1113+F1233+F1333)*a4;
C36 = 4*F1233*a3+2*(F1112+F1222+F1233)*a4;
C44 = 8*F2233*a3+4*(F1123+F2223+F2333)*a5;
C45 = 8*F1233*a3;
C46 = 8*F1233*a3;
C55 = 8*F1133*a3+4*(F1113+F1233+F1333)*a5;
C56 = 8*F1123*a3;
C66 = 8*F1122*a3+4*(F1112+F1222+F1233)*a5;
C = [C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66];