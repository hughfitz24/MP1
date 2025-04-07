%% MP1.m
% Main file implementing the user-defined functions for the solution of
% the presented complementary symmetry, push-pull power amplifier.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of Minor
% Project 1.

%% Diode Curve Fitting

clear;
clc;

V_t = 26e-3;
I_S_DIODE = 15.2e-12;

% Data Point Selection

V_data = [50e-3 100e-3 200e-3 350e-3 420e-3 480e-3 520e-3 580e-3 620e-3 680e-3 700e-3 740e-3 ...
    760e-3 780e-3 800e-3 840e-3 880e-3 920e-3 940e-3 980e-3];

I_data = [30e-12 123e-12 1.323e-9 40.4e-9 197e-9 784e-9 1.94e-6 7.56e-6 18.88e-6 73.6e-6 ...
    116e-6 284e-6 440e-6 672e-6 1.024e-3 2.240e-3 4.480e-3 8.053e-3 10.29e-3 15.52e-3];

N = length(V_data);

alpha = zeros(1, N);

for i = 1:N
    alpha(i) = V_t*log(I_data(i)/I_S_DIODE - 1);
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
hLine = plot(V_data, I_data, '-', 'LineWidth', 2);
hold on
hPoints = plot(V_data, I_data, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'r');
xlabel('Voltage (V)')
ylabel('Current (A)')
title('I-V Characteristic for IN4148 Diode, with data points highlighted')
grid on
legend([hLine, hPoints], {'I-V Characteristic', 'Chosen Data Points'}, 'Location', 'best')

figure;
plot(V, I_data, 'rx', 'DisplayName', 'Curve Fitted Data', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot(V_data, I_data, 'b-', 'DisplayName', 'Experimental Data', 'LineWidth', 1.5);
grid on;
xlabel('Voltage (V)');
ylabel('Current (A)');
title('I-V Characteristic for IN4148 Diode');
legend('Location', 'northwest');
hold off;

%% DC Operating Point

% Define knwon circuit parameters

I_S1 = 2.584e-15;
beta_F1 = 44.8;
V_A1 = 61.8682;

I_S2 = 1.871e-15;
beta_F2 = 28;
V_A2 = 123.7323;

I_S3 = 4e-12;
beta_F3 = 84;
V_A3 = 92.799;


% Define anonymous function for system of equations F

Fvec = @(x) [
    x(5) - (IS1/BF1)*(exp((x(1)-9)/VT) - 1);
    x(6) - BF1*x(5)*(1 - (x(1)-x(2))/VA1);
    x(7) - (IS2/BF2)*(exp((x(2)-x(3))/VT) - 1);
    x(8) - BF2*x(7)*(1 - x(3)/VA2);
    x(9) - (IS3/BF3)*(exp((x(4)-x(1))/VT) - 1);
    x(10) - BF3*x(9)*(1 - x(4)/VA3);
    0 - x(9) + x(7) + x(8) - x(5) + x(6);
    x(10) + x(11) + x(7) - x(9);
    (9 - x(1))/1025 - x(5) + x(11);
    x(1) - x(3) - 2*(n*VT*log(x(11)/ISD + 1) + R*x(11));
    (x(2) - x(4))/1e5 - x(9)
    ];

% Define anonymous function for Jacobian vector
Jvec = @(x) [
    -(IS1/BF1)*exp((x(1)-9)/VT),               0,                                  0,                                  0,                         1,       0,       0,       0,        0,       0,   0;
    (BF1*x(5))/VA1,                          -(BF1*x(5))/VA1,                      0,                                  0,                        -1,       1,       0,       0,        0,       0,   0;
    0,                                      -(IS2/BF2)*exp((x(2)-x(3))/VT),       (IS2/BF2)*exp((x(2)-x(3))/VT),       0,                         0,       0,       1,       0,        0,       0,   0;
    0,                                       0,                                   (BF2*x(7))/VA2,                      0,                         0,       0,  -BF2*(1 - x(3)/VA2), 1,        0,       0,   0;
    (IS3/BF3)*exp((x(4)-x(1))/VT),           0,                                   0,                          -(IS3/BF3)*exp((x(4)-x(1))/VT),  0,       0,       0,       0,        1,       0,   0;
    0,                                       0,                                   0,                          (BF3*x(9))/VA3,                  0,       0,       0,       0,  -BF3*(1 - x(4)/VA3), 1, 0;
    0,                                       0,                                   0,                                  0,                        -1,       1,       1,       1,       -1,       0,   0;
    0,                                       0,                                   0,                                  0,                         0,       0,      -1,       0,        1,       1,   1;
    -1/1025,                                  0,                                   0,                                  0,                        -1,       0,       0,       0,        0,       0,   1;
    1,                                       0,                                  -1,                                  0,                         0,       0,       0,       0,        0,       0,  -2*n*VT*(1/(x(11)+ISD) + R);
    0,                                   1/1e5,                                   0,                            -1/1e5,                         0,       0,       0,       0,       -1,       0,   0
    ];