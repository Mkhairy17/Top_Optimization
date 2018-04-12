function [C] = Stiffness_Matrix(a1,a2,a3,a4,a5,F1111,F1122,F1133,F2222,F2233,F3333,F2333,F2223,F1123,F1113,F1233,F1333,F1112,F1222,F1223)
C1111 = 2*a1+2*a2+2*F1111*a3+2*(F1111+F1122+F1133)*a4+2*(F1111+F1122+F1133)*a5; 
C1122 = 2*a1+2*F1122*a3+(F1111+F1122+F1133)*a4+(F1122+F2222+F2233)*a4;
C1133 = 2*a1+2*F1133*a3+(F1111+F1122+F1133)*a4+(F1133+F2233+F3333)*a4;
C1123 = 4*F1123*a3+2*(F1123+F2223+F2333)*a4;
C1113 = 4*F1113*a3+2*(F1113+F1233+F1333)*a4;
C1112 = 4*F1112*a3+2*(F1112+F1222+F1233)*a4;
C2222 = 2*a1+2*a2+2*F2222*a3+2*(F1122+F2222+F2233)*a4+2*(F1122+F2222+F2233)*a5;
C2233 = 2*a1+2*F2233*a3+(F1122+F2222+F2233)*a4+(F1133+F2233+F3333)*a4;
C2223 = 4*F2223*a3+2*(F1123+F2223+F2333)*a4;
C2213 = 4*F1223*a3+2*(F1113+F1233+F1333)*a4;
C2212 = 4*F1222*a3+2*(F1112+F1222+F1233)*a4;
C3333 = 2*a1+2*a2+2*F3333*a3+2*(F1133+F2233+F3333)*a4+2*(F1133+F2233+F3333)*a5;
C3323 = 4*F2333*a3+2*(F1123+F2223+F2333)*a4;
C3313 = 4*F1333*a3+2*(F1113+F1233+F1333)*a4;
C3312 = 4*F1233*a3+2*(F1112+F1222+F1233)*a4;
C2323 = 8*F2233*a3+4*(F1123+F2223+F2333)*a5;
C2313 = 8*F1233*a3;
C2312 = 8*F1233*a3;
C1313 = 8*F1133*a3+4*(F1113+F1233+F1333)*a5;
C1312 = 8*F1123*a3;
C1212 = 8*F1122*a3+4*(F1112+F1222+F1233)*a5;
C = [C1111 C1122 C1133 C1123 C1113 C1112 ;C1122 C2222 C2233 C2223 C2213 C2212;C1133 C2233 C3333 C3323 C3313 C3312;C1123 C2223 C3323 C2323 C2313 C2312;C1113 C2213 C3313 C2313 C1313 C1312;C1112 C2212 C3312 C2312 C1312 C1212];

