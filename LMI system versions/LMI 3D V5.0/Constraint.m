function F1 = Constraint (V1 ,V2 ,V3 ,V4 ,V5 ,V6 ,V7 ,V8 ,V9, V10 ,V11 ,V12, V13, V14)
F1=[3*(56+120*V1-20*V10+3*V14-40*V5+35*V6),56+3*V14-40*V5-105*V6,4*(14+15*V1+15*V10-3*V14+5*V5),60*(V13+4*V4-7*V9),60*(-3*V12+12*V3+7*V8),30*(2*V11+24*V2-7*V7);... 
56+3*V14-40*V5-105*V6,3*(56-120*V1+20*V10+3*V14-40*V5+35*V6),4*(14-15*V1-15*V10-3*V14+5*V5),60*(3*V13+12*V4+7*V9),60*(-V12+4*V3-7*V8),30*(2*V11+24*V2+7*V7);...
4*(14+15*V1+15*V10-3*V14+5*V5),4*(14-15*V1-15*V10-3*V14+5*V5),24*(7+V14+10*V5),240*(-V13+3*V4),240*(V12+3*V3),120*(-V11+2*V2);...
60*(V13+4*V4-7*V9),60*(3*V13+12*V4+7*V9),240*(-V13+3*V4),16*(14-15*V1-15*V10-3*V14+5*V5),240*(-V11+2*V2),120*(-V12+4*V3-7*V8);...
60*(-3*V12+12*V3+7*V8),60*(-V12+4*V3-7*V8),240*(V12+3*V3),240*(-V11+2*V2),16*(14+15*V1+15*V10-3*V14+5*V5),120*(V13+4*V4-7*V9);...
30*(2*V11+24*V2-7*V7),30*(2*V11+24*V2+7*V7),120*(-V11+2*V2),120*(-V12+4*V3-7*V8),120*(V13+4*V4-7*V9),4*(56+3*V14-40*V5-105*V6)]