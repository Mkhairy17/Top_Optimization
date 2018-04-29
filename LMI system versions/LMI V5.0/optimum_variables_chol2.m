function xval_opt = optimum_variables_chol2 (U1,U2,U3,U4,U5,L)

cvx_begin sdp quiet
    variables F11 F12 F13 F22 F23 F33 V1 V2 V3 V4
    F = [F11 F12 F13;F12 F22 F23;F13 F23 F33];
%     V = [V1;V2;V3;V4];
    A11 =  U1 + U2.*V1 + U3.*V3;
    A12 =  U4 - U3.*V3;  
    A16 =  U2*V2/2 + U3*V4;
    A26 =  U2*V2/2 - U3*V4;
    A22 =  U1 - U2.*V1 + U3.*V3;
    A66 =  U5 - U3.*V3;
    A = [A11 A12 A16;A12 A22 A26;A16 A26 A66];
    minimize( F11+F22+F33 );
   [A L';L F] >=0;
    [3+4*V1+V3,4*V2+2*V4,1-V3;4*V2+2*V4,4-4*V3,4*V2-2*V4;1-V3,4*V2-2*V4,3-4*V1+V3] >= 0;
cvx_end
xval_opt = [V1 V2 V3 V4];
