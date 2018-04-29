function xval_opt = optimum_variables2(a1,a2,a3,a4,a5,L)
cvx_begin sdp quiet
    variables FF11 FF12 FF13 FF14 FF15 FF16 FF22 FF23 FF24 FF25 FF26 FF33 FF34 FF35 FF36 FF44 FF45 FF46 FF55 FF56 FF66 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
    F = [FF11 FF12 FF13 FF14 FF15 FF16;...
         FF12 FF22 FF23 FF24 FF25 FF26;...
         FF13 FF23 FF33 FF34 FF35 FF36;...
         FF14 FF24 FF34 FF44 FF45 FF46;...
         FF15 FF25 FF35 FF45 FF55 FF56;...
         FF16 FF26 FF36 FF46 FF56 FF66];
    F1111 = 1/280 *(56 + 120 *V1 - 20 *V10 + 3 *V14 - 40 *V5 + 35 *V6);
    F1122 = 1/840 * (56 + 3 *V14 - 40* V5 - 105 *V6);
    F1133 = 1/210 * (14 + 15 *V1 + 15 *V10 - 3 *V14 + 5* V5);
    F1123 = 1/28 *(V13 + 4 *V4 - 7 *V9);
    F1113 = 1/28 *(-3 * V12 + 12* V3 + 7 *V8);
    F1112 = 1/56 *(2* V11 + 24 * V2 - 7 * V7);
    F2222 = 1/280 *(56 - 120 *V1 + 20 *V10 + 3 *V14 - 40* V5 + 35 *V6);
    F2233 = 1/210 *(14 - 15 *V1 - 15 * V10 - 3 * V14 + 5 * V5);
    F2223 = 1/28 *(3* V13 + 12 *V4 + 7 * V9);
    F1223 = 1/28 *(-V12 + 4*  V3 - 7 * V8);
    F1222 = 1/56 *(2 * V11 + 24* V2 + 7* V7);
    F3333 = 1/35 *(7 + V14 + 10 *V5);
    F2333 = 1/7 *(-V13 + 3 *V4);
    F1333 = 1/7 *(V12 + 3 *V3);
    F1233 = 1/14 *(-V11 + 2 * V2);
   %%
    C11 = 2*a1+ 2*a2+2*F1111*a3+ 2*(F1111+F1122+F1133)*a4+ 2*(F1111+F1122+F1133)*a5;
    C12 = 2*a1+ 2*F1122*a3+ (F1111+F1122+F1133)*a4+ (F1122+F2222+F2233)*a4;
    C13 = 2*a1+2*F1133*a3+(F1111+F1122+F1133)*a4+(F1133+F2233+F3333)*a4;
    C14 = 4*F1123*a3+ 2*(F1123+F2223+F2333)*a4;
    C15 = 4*F1113*a3+2*(F1113+F1223+F1333)*a4; 
    C16 = 4*F1112*a3+ 2*(F1112+F1222+F1233)*a4;
    C22 = 2*a1+ 2*a2+ 2*F2222*a3+ 2*(F1122+F2222+F2233)*a4+ 2*(F1122+F2222+F2233)*a5;
    C23 = 2*a1+ 2*F2233*a3+ (F1122+F2222+F2233)*a4+ (F1133+F2233+F3333)*a4;
    C24 = 4*F2223*a3+ 2*(F1123+F2223+F2333)*a4;
    C25 = 4*F1223*a3+2*(F1113+F1223+F1333)*a4;
    C26 = 4*F1222*a3+ 2*(F1112+F1222+F1233)*a4;
    C33 = 2*a1+2*a2+2*F3333*a3+2*(F1133+F2233+F3333)*a4+2*(F1133+F2233+F3333)*a5;
    C34 = 4*F2333*a3+ 2*(F1123+F2223+F2333)*a4;
    C35 = 4*F1333*a3+2*(F1113+F1223+F1333)*a4;
    C36 = 4*F1233*a3+ 2*(F1112+F1222+F1233)*a4;
    C44 = 8*F2233*a3+ 4*(F1123+F2223+F2333)*a5;
    C45 = 8*F1233*a3;
    C46 = 8*F1223*a3;
    C55 = 8*F1133*a3+4*(F1113+F1223+F1333)*a5;
    C56 = 8*F1123*a3;
    C66 = 8*F1122*a3+ 4*(F1112+F1222+F1233)*a5;
    C = [C11 C12 C13 C14 C15 C16;C12 C22 C23 C24 C25 C26;C13 C23 C33 C34 C35 C36;C14 C24 C34 C44 C45 C46;C15 C25 C35 C45 C55 C56;C16 C26 C36 C46 C56 C66];  
 %%
    minimize( trace(F) )
    [C L';L F] >=0
    [3*(56+120*V1-20*V10+3*V14-40*V5+35*V6),56+3*V14-40*V5-105*V6,4*(14+15*V1+15*V10-3*V14+5*V5),60*(V13+4*V4-7*V9),60*(-3*V12+12*V3+7*V8),30*(2*V11+24*V2-7*V7);... 
    56+3*V14-40*V5-105*V6,3*(56-120*V1+20*V10+3*V14-40*V5+35*V6),4*(14-15*V1-15*V10-3*V14+5*V5),60*(3*V13+12*V4+7*V9),60*(-V12+4*V3-7*V8),30*(2*V11+24*V2+7*V7);...
    4*(14+15*V1+15*V10-3*V14+5*V5),4*(14-15*V1-15*V10-3*V14+5*V5),24*(7+V14+10*V5),240*(-V13+3*V4),240*(V12+3*V3),120*(-V11+2*V2);...
    60*(V13+4*V4-7*V9),60*(3*V13+12*V4+7*V9),240*(-V13+3*V4),16*(14-15*V1-15*V10-3*V14+5*V5),240*(-V11+2*V2),120*(-V12+4*V3-7*V8);...
    60*(-3*V12+12*V3+7*V8),60*(-V12+4*V3-7*V8),240*(V12+3*V3),240*(-V11+2*V2),16*(14+15*V1+15*V10-3*V14+5*V5),120*(V13+4*V4-7*V9);...
    30*(2*V11+24*V2-7*V7),30*(2*V11+24*V2+7*V7),120*(-V11+2*V2),120*(-V12+4*V3-7*V8),120*(V13+4*V4-7*V9),4*(56+3*V14-40*V5-105*V6)] >= 0
cvx_end
xval_opt=[V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14];