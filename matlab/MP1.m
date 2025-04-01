%% MP1.m
% Main file implementing the user-defined functions for the solution of
% the presented complementary symmetry, push-pull power amplifier.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of Minor
% Project 1.

%% Diode Curve Fitting

clear;
clc;

V_t = 26e-3;
Is = 15.2e-12;

% Data Point Selection

V_data = [50e-3 100e-3 200e-3 350e-3 420e-3 480e-3 520e-3 580e-3 620e-3 680e-3 700e-3 740e-3 ...
    760e-3 780e-3 800e-3 840e-3 880e-3 920e-3 940e-3 980e-3];

I_data = [30e-12 123e-12 1.323e-9 40.4e-9 197e-9 784e-9 1.94e-6 7.56e-6 18.88e-6 73.6e-6 ...
    116e-6 284e-6 440e-6 672e-6 1.024e-3 2.240e-3 4.480e-3 8.053e-3 10.29e-3 15.52e-3];

N = length(V_data);

alpha = zeros(1, N);

for i = 1:N
    alpha(i) = V_t*log(I_data(i)/Is - 1);
end

beta = I_data;

A = [sum(alpha.^2) sum(alpha.*beta);
     sum(alpha.*beta) sum(beta.^2)];

b = [-sum(V_data.*(-alpha));
-sum(V_data.*(-beta))];

x = GaussianElimination(A,b);

n = x(1);
R = x(2);

MSE = (1/N) * sum((V_data - (alpha .* n + beta .* R)).^2);

fprintf("n = %.5f\n", n)
fprintf("R = %.5f\n", R)
fprintf("MSE = %.10f", MSE)

V = alpha*n + beta*R;

figure
plot(V_data, I_data, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'r')
xlabel('Voltage (V)')
ylabel('Current (A)')
title('I-V Characteristic for IN4148 Diode, with data points highlighted')
grid on

figure
plot(V, I_data, 'rx')
hold on;
plot(V_data, I_data, 'b-')