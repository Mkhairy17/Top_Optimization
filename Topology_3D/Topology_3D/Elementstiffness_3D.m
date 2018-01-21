function [KE] = Elementstiffness_3D(a,b,c)

syms x y z 
E = 1;
v = 0.3;
c11 = 1/ E;
c12 = -v/E;
c13 = -v/E;
c14 = 0;
c15 = 0;
c16 = 0;
c22 = 1/E;
c23 = -v/E;
c24 = 0;
c25 = 0;
c26 = 0;
c33 = 1/E;
c34 = 0;
c35 = 0;
c36 = 0;
c44 = 2*(1+v)/E;
c45 = 0;
c46 = 0;
c55 = 2*(1+v)/E;
c56 = 0;
c66 = 2*(1+v)/E;

N = 1/8*[(1-2*x/a)*(1-2*y/b)*(1-2*z/c);(1+2*x/a)*(1-2*y/b)*(1-2*z/c);(1+2*x/a)*(1+2*y/b)*(1-2*z/c);(1-2*x/a)*(1+2*y/b)*(1-2*z/c);(1-2*x/a)*(1-2*y/b)*(1+2*z/c);(1+2*x/a)*(1-2*y/b)*(1+2*z/c);(1+2*x/a)*(1+2*y/b)*(1+2*z/c);(1-2*x/a)*(1+2*y/b)*(1+2*z/c)];
B11 = diff(N(1),x);
B12 = diff(N(1),y);
B13 = diff(N(1),z);
B21 = diff(N(2),x);
B22 = diff(N(2),y);
B23 = diff(N(2),z);
B31 = diff(N(3),x);
B32 = diff(N(3),y);
B33 = diff(N(3),z);
B41 = diff(N(4),x);
B42 = diff(N(4),y);
B43 = diff(N(4),z);
B51 = diff(N(5),x);
B52 = diff(N(5),y);
B53 = diff(N(5),z);
B61 = diff(N(6),x);
B62 = diff(N(6),y);
B63 = diff(N(6),z);
B71 = diff(N(7),x);
B72 = diff(N(7),y);
B73 = diff(N(7),z);
B81 = diff(N(8),x);
B82 = diff(N(8),y);
B83 = diff(N(8),z);

B = [B11 0 0 B21 0 0 B31 0 0 B41 0 0 B51 0 0 B61 0 0 B71 0 0 B81 0 0;...
    0 B12 0 0 B22 0 0 B32 0 0 B42 0 0 B52 0 0 B62 0 0 B72 0 0 B82 0;...
    0 0 B13 0 0 B23 0 0 B33 0 0 B43 0 0 B53 0 0 B63 0 0 B73 0 0 B83;...
    0 B13 B12 0 B23 B22 0 B33 B32 0 B43 B42 0 B53 B52 0 B63 B62 0 B73 B72 0 B83 B82;...
    B13 0 B11 B23 0 B21 B33 0 B31 B43 0 B41 B53 0 B51 B63 0 B61 B73 0 B71 B83 0 B81;...
    B12 B11 0 B22 B21 0 B32 B31 0 B42 B41 0 B52 B51 0 B62 B61 0 B72 B71 0 B82 B81 0
    ];

s_6 = [c11 c12 c13 c14 c15 c16;c12 c22 c23 c24 c25 c26;c13 c23 c33 c34 c35 c36;c14 c24 c34 c44 c45 c46;c15 c25 c35 c45 c55 c56;c16 c26 c36 c46 c56 c66];
C = inv(s_6);

ke1 = int(B'*C*B,x,-a/2,a/2);
ke2 = int(ke1,y,-b/2,b/2);
KE = double(int(ke2,z,-c/2,c/2));