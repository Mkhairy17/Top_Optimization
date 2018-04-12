function xval_opt = optimum_variables (U1,U2,U3,U4,U5,vv1,vv2,vv3)
c = [1 1 1 0 0 0 0]';
setlmis([]);
alpha1=lmivar(2,[1,1]); 
alpha2=lmivar(2,[1,1]);
alpha3=lmivar(2, [1,1]);
V1 = lmivar(2,[1,1]);
V2 = lmivar(2,[1,1]);
V3 = lmivar(2,[1,1]);
V4 = lmivar(2,[1,1]);
%% LMI 1 (A V1;V1' alfa1) >= 0
lmiterm([-1 4 1 0], vv1(1));   %V1'
lmiterm([-1 4 2 0], vv1(2));   %V1'
lmiterm([-1 4 3 0], vv1(3));   %V1'
lmiterm([-1 4 4 alpha1],1,1);   %alfa1
% lmiterm([-1 1 4 0],vv1(1));   %V1'
% lmiterm([-1 2 4 0],vv1(2));   %V1'
% lmiterm([-1 3 4 0],vv1(3));   %V1'
% A11 =  U1 + U2.*V1 + U3.*V3;
lmiterm([-1 1 1 0], U1);  %U1
lmiterm([-1 1 1 V1], U2 , 1);  %U2.*V1
lmiterm([-1 1 1 V3],U3,1);  %U3.*V3
% A12 =  U4 - U3.*V3;
lmiterm([-1 1 2 0],U4);  %U4
lmiterm([-1 1 2 V3],U3,-1);  %U3.*V3

%A16 =  U2.*V2/2 + U3.*V4;
lmiterm([-1 1 3 V2],U2/2,1);  % U2.*V2/2
lmiterm([-1 1 3 V4],U3,1);  % U3.*V4

%A26 =  U2.*V2/2 - U3.*V4;
lmiterm([-1 2 3 V2],U2/2,1);  % U2.*V2/2
lmiterm([-1 2 3 V4],U3,-1);  % U3.*V4

%A22 =  U1 - U2.*V1 + U3.*V3;
lmiterm([-1 2 2 0],U1);  %U1
lmiterm([-1 2 2 V1],U2,-1);  % U2.*V1
lmiterm([-1 2 2 V3],U3, 1);  % U3.*V3

%A66 =  U5 - U3.*V3;

lmiterm([-1 3 3 0], U5);  %U5
lmiterm([-1 3 3 V3],U3,-1);  % U3.*V3
%==========================================

%% LMI 2 (A V2;V2' alfa2) >= 0

lmiterm([-2 4 1 0], vv2(1));   %V1'
lmiterm([-2 4 2 0], vv2(2));   %V1'
lmiterm([-2 4 3 0], vv2(3));   %V1'
lmiterm([-2 4 4 alpha2],1,1);   %alfa1
% lmiterm([-2 1 4 0],vv2(1));   %V1'
% lmiterm([-2 2 4 0],vv2(2));   %V1'
% lmiterm([-2 3 4 0],vv2(3));   %V1'

% A11 =  U1 + U2.*V1 + U3.*V3;
lmiterm([-2 1 1 0],U1);  %U1
lmiterm([-2 1 1 V1],U2, 1);  %U2.*V1
lmiterm([-2 1 1 V3],U3, 1);  %U3.*V3
% A12 =  U4 - U3.*V3;
lmiterm([-2 1 2 0], U4);  %U4
lmiterm([-2 1 2 V3],U3 ,-1);  %U3.*V3

%A16 =  U2.*V2/2 + U3.*V4;
lmiterm([-2 1 3 V2], U2/2 , 1);  % U2.*V2/2
lmiterm([-2 1 3 V4], U3 , 1);  % U3.*V4

%A26 =  U2.*V2/2 - U3.*V4;
lmiterm([-2 2 3 V2],U2/2 , 1);  % U2.*V2/2
lmiterm([-2 2 3 V4],U3 ,-1);  % U3.*V4

%A22 =  U1 - U2.*V1 + U3.*V3;
lmiterm([-2 2 2 0],U1);  %U1
lmiterm([-2 2 2 V1],U2 ,-1);  % U2.*V1
lmiterm([-2 2 2 V3], U3 , 1);  % U3.*V3

%A66 =  U5 - U3.*V3;

lmiterm([-2 3 3 0], U5);  %U5
lmiterm([-2 3 3 V3], U3 ,-1);  % U3.*V3
%============================================

%% LMI 3 (A V3;V3' alfa3) >= 0

lmiterm([-3 4 1 0], vv3(1));   %V1'
lmiterm([-3 4 2 0], vv3(2));   %V1'
lmiterm([-3 4 3 0], vv3(3));   %V1'
lmiterm([-3 4 4 alpha3],1,1);   %alfa1
% lmiterm([-3 1 4 0],vv3(1));   %V1'
% lmiterm([-3 2 4 0],vv3(2));   %V1'
% lmiterm([-3 3 4 0],vv3(3));   %V1'

% A11 =  U1 + U2.*V1 + U3.*V3;
lmiterm([-3 1 1 0], U1);  %U1
lmiterm([-3 1 1 V1], U2 , 1);  %U2.*V1
lmiterm([-3 1 1 V3], U3 , 1);  %U3.*V3
% A12 =  U4 - U3.*V3;
lmiterm([-3 1 2 0], U4);  %U4
lmiterm([-3 1 2 V3], U3 ,-1);  %U3.*V3

%A16 =  U2.*V2/2 + U3.*V4;
lmiterm([-3 1 3 V2], U2/2, 1);  % U2.*V2/2
lmiterm([-3 1 3 V4], U3, 1);  % U3.*V4

%A26 =  U2.*V2/2 - U3.*V4;
lmiterm([-3 2 3 V2], U2/2, 1);  % U2.*V2/2
lmiterm([-3 2 3 V4], U3,-1);  % U3.*V4

%A22 =  U1 - U2.*V1 + U3.*V3;
lmiterm([-3 2 2 0], U1);  %U1
lmiterm([-3 2 2 V1], U2, -1);  % U2.*V1
lmiterm([-3 2 2 V3], U3, 1);  % U3.*V3

%A66 =  U5 - U3.*V3;

lmiterm([-3 3 3 0], U5);  %U5
lmiterm([-3 3 3 V3], U3,-1);  % U3.*V3

%==========================================

% LMI 4 F >= 0
lmiterm([-4 1 1 0], 3);   %CONSTANT = 3
lmiterm([-4 1 1 V1], 4,1);   %4V1
lmiterm([-4 1 1 V3], 1,1);   % V3

lmiterm([-4 1 2 V2], 4,1);   % 4V2 
lmiterm([-4 1 2 V4], 2,1);   % 2V4


lmiterm([-4 1 3 0], 1);   % CONSTANT = 1 
lmiterm([-4 1 3 V3], -1,1);   % -1 * V3


lmiterm([-4 2 2 0], 4);   % CONSTANT = 4 
lmiterm([-4 2 2 V3], 4,-1);   % -4 * V3

lmiterm([-4 2 3 V2], 4,1);   % 4 V2 
lmiterm([-4 2 3 V4], 2,-1);   % -2 V4


lmiterm([-4 3 3 0],3);   % CONSTANT = 3 
lmiterm([-4 3 3 V1],-4,1);   % - 4 V2 
lmiterm([-4 3 3 V3],1,1);   %  V3


lmiterm([-4 2 1 V2], 4,1);   % 4V2 
lmiterm([-4 2 1 V4], 2,1);   % 2V4

lmiterm([-4 3 1 0], 1);   % CONSTANT = 1 
lmiterm([-4 3 1 V3], -1,1);   % -1 * V3


lmiterm([-4 3 2 V2], 4,1);   % 4 V2 
lmiterm([-4 3 2 V4], -2,1);   % -2 V4



%=========================================

%Create the LMI system
lmisys = getlmis;
options = [1e-3, 100, -1, 5, 1];
[copt, xopt] = mincx (lmisys, c , options);
xval_opt = xopt;
evals = evallmi(lmisys,xopt);
% [lhs1,rhs1] = showlmi(evals,1) 
% [lhs2,rhs2] = showlmi(evals,2)
% [lhs2,rhs2] = showlmi(evals,3)
% [lhs2,rhs2] = showlmi(evals,4)

