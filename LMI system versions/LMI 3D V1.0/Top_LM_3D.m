clc
clear all
%% User inputs
Lx = input('length of x=');
Ly = input('length of y=');
Lz = input('length of z=');
nx = input('nx=');
ny = input('ny=');
nz = input('nz=');
ne = nx*ny*nz;
a = Lx/nx; %elemeny width
b = Ly/ny; %element length
c = Lz/nz; %element height
%%
%initial values
xval_old = zeros(15*Ne,1);
xval_opt = zeros(7,Ne);
% A Matrix should be defined here
%%
% Displacement Solver



