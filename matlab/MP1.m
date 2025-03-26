%% MP1.m
% Main file implementing the user-defined functions for the solution of
% the presented complementary symmetry, push-pull power amplifier.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of Minor 
% Project 1. 

%% Diode Curve Fitting

% Initial Guesses

n_init = 1;
R_init = 8.522425;

% Data Point Selection

V_data = [10e-3 50e-3 200e-3 350e-3 420e-3 480e-3 520e-3 580e-3 620e-3 680e-3 700e-3 740e-3 ...
    760e-3 780e-3 800e-3 840e-3 880e-3 920e-3 940e-3 980e-3];

I_data = [3.73e-12 30e-12 1.323e-9 40.4e-9 197e-9 784e-9 1.94e-6 7.56e-6 18.88e-6 73.6e-6 ...
     116e-6 284e-6 440e-6 672e-6 1.024e-3 2.240e-3 4.480e-3 8.053e-3 10.29e-3 15.52e-3];

J = @(x) [Vt, Id];
V = @(Id, n, R) ln(Id/Is+1) + n*Vt + R*Id;

