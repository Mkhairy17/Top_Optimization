function xval_opt = optimum_variables_chol (U1,U2,U3,U4,U5,L)
c = [1 0 0 1 0 1 0 0 0 0]';
setlmis([]);
F11=lmivar(2,[1,1]); 
F12=lmivar(2,[1,1]);
F13=lmivar(2, [1,1]);
F22=lmivar(2,[1,1]);
F23=lmivar(2, [1,1]);
F33=lmivar(2, [1,1]);
V1 = lmivar(2,[1,1]);
V2 = lmivar(2,[1,1]);
V3 = lmivar(2,[1,1]);
V4 = lmivar(2,[1,1]);
%% LMI 1 (A V1;V1' alfa1) >= 0
lmiterm([-1 4 1 0], L(1,1));   %V1'
lmiterm([-1 4 2 0], L(1,2));   %V1'
lmiterm([-1 4 3 0], L(1,3));   %V1'
lmiterm([-1 4 4 F11],1,1);   %V1'
lmiterm([-1 4 5 F12],1,1);   %V1'
lmiterm([-1 4 6 F13],1,1);   %V1'

lmiterm([-1 5 1 0], L(2,1));   %V1'
lmiterm([-1 5 2 0], L(2,2));   %V1'
lmiterm([-1 5 3 0], L(2,3));   %V1'
lmiterm([-1 5 4 F12],1,1);   %V1'
lmiterm([-1 5 5 F22],1,1);   %V1'
lmiterm([-1 5 6 F23],1,1);   %V1'

lmiterm([-1 6 1 0], L(3,1));   %V1'
lmiterm([-1 6 2 0], L(3,2));   %V1'
lmiterm([-1 6 3 0], L(3,3));   %V1'
lmiterm([-1 6 4 F13],1,1);   %V1'
lmiterm([-1 6 5 F23],1,1);   %V1'
lmiterm([-1 6 6 F33],1,1);   %V1'

%%
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

% LMI 4 F >= 0
lmiterm([-2 1 1 0], 3);   %CONSTANT = 3
lmiterm([-2 1 1 V1], 4,1);   %4V1
lmiterm([-2 1 1 V3], 1,1);   % V3

lmiterm([-2 1 2 V2], 4,1);   % 4V2 
lmiterm([-2 1 2 V4], 2,1);   % 2V4


lmiterm([-2 1 3 0], 1);   % CONSTANT = 1 
lmiterm([-2 1 3 V3], -1,1);   % -1 * V3


lmiterm([-2 2 2 0], 4);   % CONSTANT = 4 
lmiterm([-2 2 2 V3], 4,-1);   % -4 * V3

lmiterm([-2 2 3 V2], 4,1);   % 4 V2 
lmiterm([-2 2 3 V4], 2,-1);   % -2 V4


lmiterm([-2 3 3 0],3);   % CONSTANT = 3 
lmiterm([-2 3 3 V1],-4,1);   % - 4 V2 
lmiterm([-2 3 3 V3],1,1);   %  V3


% lmiterm([-2 2 1 V2], 4,1);   % 4V2 
% lmiterm([-2 2 1 V4], 2,1);   % 2V4
% 
% lmiterm([-2 3 1 0], 1);   % CONSTANT = 1 
% lmiterm([-2 3 1 V3], -1,1);   % -1 * V3


% lmiterm([-2 3 2 V2], 4,1);   % 4 V2 
% lmiterm([-2 3 2 V4], -2,1);   % -2 V4



%=========================================
%Create the LMI system
lmisys = getlmis;
options = [1e-9, 100, -1, 5, 1];
[copt, xopt] = mincx (lmisys, c , options);
xval_opt = xopt;
evals = evallmi(lmisys,xopt);
[lhs1,rhs1] = showlmi(evals,1) ;



