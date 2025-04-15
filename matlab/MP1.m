%% MP1.m
% Main file implementing the user-defined functions for the solution of
% the presented complementary symmetry, push-pull power amplifier.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of Minor
% Project 1.

clear;
clc;

V_t = 26e-3;

%% Diode Curve Fitting

I_S_DIODE = 15.2e-12; % A

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

x_gauss= GaussianElimination(A,b);

n_gauss = x_gauss(1);
R_gauss= x_gauss(2);

MSE_gauss = (1/N) * sum((V_data - (alpha .* n_gauss + beta .* R_gauss)).^2);

fprintf("n = %.5f\n", n_gauss)
fprintf("R = %.5f\n", R_gauss)
fprintf("MSE = %.10f\n", MSE_gauss)

% Solve with Gauss Seidel
x_GS = GaussSeidel(A, b, 1e-9, 1000, [0; 0]);
n_GS = x_GS(1);
R_GS = x_GS(2);

MSE_GS = (1/N) * sum((V_data - (alpha .* n_GS + beta .* R_GS)).^2);

fprintf("n = %.5f\n", n_GS)
fprintf("R = %.5f\n", R_GS)
fprintf("MSE = %.10f\n", MSE_GS)

if (MSE_gauss < MSE_GS)
    fprintf("Gauss Elimination found more accurate solution. \n")
    n = n_gauss;
    R = R_gauss;
    MSE = MSE_gauss;
else
    fprintf("Gauss-Seidel Iteration found more accurate solution. \n")
    n = n_GS;
    R = R_GS;
    MSE = MSE_GS;
end

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

% Define known circuit parameters

I_S_DIODE = 15.2e-9; % mA, scaled for convergence
V_t = 26e-3; % V

I_S1 = 2.6e-12; %mA
beta_F1 = 44.8;
V_A1 = 61.87;

I_S2 = 1.87e-12; %mA
beta_F2 = 28;
V_A2 = 123.73;

I_S3 = 4e-9; %mA
beta_F3 = 84;
V_A3 = 92.8;

R = R/1000; % Convert to kOhms

% Define anonymous function for system of equations F

% NOTE: F takes a vector x as input, where:
% x(1) = VB1
% x(2) = VE1
% x(3) = VB2
% x(4) = VB3
% x(5) = IB1
% x(6) = IC1
% x(7) = IB2
% x(8) = IC2
% x(9) = IB3
% x(10) = IC3
% x(11) = ID

F = @(x) [
    x(5) - (I_S1/beta_F1)*(exp((x(1) - x(2))/V_t) - 1);                                 % F1
    x(6) - (beta_F1 * x(5))*(1 - (x(1) - 9)/V_A1);                                           % F2
    x(7) - (I_S2/beta_F2)*(exp((x(2) - x(3))/V_t) - 1);                                 % F3
    x(8) - beta_F2*x(7)*(1 -(0 - x(3))/V_A2);                                               % F4
    x(9) - (I_S3/beta_F3)*(exp((x(4))/V_t) - 1);                                 % F5
    x(10) - beta_F3*x(9)*(1 - (x(4)-x(3))/V_A3);                                              % F6
    x(5) + x(6) - x(7) - x(8) - x(9);                                             % F7
    x(11) + x(7) - x(10);                                          % F8
    (9 - x(1))/1.025 - x(5) - x(11);                                        % F9
    x(1) - x(3) - 2*(n*V_t*log(x(11)/I_S_DIODE + 1) - R*x(11));                                                      % F10
    (x(2) - x(4))/100 - x(9);                   % F11
    ];

% Define anonymous function for Jacobian vector

J = @(x) [
    % F1
    [-I_S1/(beta_F1*V_t)*exp((x(1)-x(2))/V_t),  I_S1/(beta_F1*V_t)*exp((x(1)-x(2))/V_t), 0, 0, 1, 0, 0, 0, 0, 0, 0];
    % F2
    [beta_F1*x(5)/V_A1, 0, 0, 0, -beta_F1 + beta_F1*(x(1)-9)/V_A1, 1, 0, 0, 0, 0, 0];
    % F3
    [0, -I_S2/(beta_F2*V_t)*exp((x(2) - x(3))/V_t), I_S2/(beta_F2*V_t)*exp((x(2) - x(3))/V_t), 0, 0, 0, 1, 0, 0, 0, 0];
    % F4
    [0, 0, -beta_F2*x(7)/V_A2, 0, 0, 0, -beta_F2*(1 + x(3)/V_A2), 1, 0, 0, 0];
    % F5
    [0, 0, 0, -I_S3/(beta_F3*V_t)*exp((x(4))/V_t), 0, 0, 0, 0, 1, 0, 0];
    % F6
    [0, 0, -beta_F3*x(9)/V_A3, beta_F3*x(9)/V_A3, 0, 0, 0, 0, -beta_F3*(1 - (x(4)-x(3))/V_A3), 1, 0];
    % F7
    [0, 0, 0, 0, 1, 1, -1, -1, -1, 0, 0];
    % F8
    [0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 1];
    % F9
    [-1/1.025, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1];
    % F10
    [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -2*(n*V_t/(x(11)+I_S_DIODE) - R)];
    % F11
    [0, 1/100, 0, -1/100, 0, 0, 0, 0, -1, 0, 0];
    ];



% Initial guess for the solution vector
x0 = [
    5.5; % VB1
    4.8; % VE1
    4.1; % VB2
    0.6; % VB3
    2.7; % IB1 (mA)
    120.96; % IC1 (mA)
    3.9; % IB2 (mA)
    109.2 % IC2 (mA)
    35e-3; % IB3 (mA)
    2.94; % IC3 (mA)
    1.15; % ID (mA)
    ];

% Tolerance and maximum iterations
tol = 1e-9;
maxIter = 20;

fprintf('Starting Newton-Raphson method...\n');

% Solve system of equations using Newton-Raphson
[gaussianSolution, gaussianIterations] = NewtonRaphson(F, J, x0, tol, maxIter);

% Display results
fprintf('--------------------------------------------------------------\n');
fprintf('Solution:\n');
fprintf('Iterations: %d\n', gaussianIterations);
fprintf('VB1 = %.4f V\n', gaussianSolution(1));
fprintf('VE1 = %.4f V\n', gaussianSolution(2));
fprintf('VB2 = %.4f V\n', gaussianSolution(3));
fprintf('VB3 = %.4f V\n', gaussianSolution(4));

fprintf('IB1 = %.9f mA\n', gaussianSolution(5));
fprintf('IC1 = %.9f mA\n', gaussianSolution(6));
fprintf('IB2 = %.9f mA\n', gaussianSolution(7));
fprintf('IC2 = %.9f mA\n', gaussianSolution(8));
fprintf('IB3 = %.9f mA\n', gaussianSolution(9));
fprintf('IC3 = %.9f mA\n', gaussianSolution(10));
fprintf('ID = %.9f mA\n', gaussianSolution(11));   
fprintf('Newton-Raphson with Gaussian Elimination completed.\n');   
fprintf('--------------------------------------------------------------\n');

%% AC Analysis

% Vector of unknown voltages and currents takes the format:

% [vb1, ve1, vb2, ib1, ic1, ib2, ic2, ib3, ic3, id]

% Small signal parameters

gm1 = 2.765;
ro1 = 860.72;
rp1 = 17.03;

gm2 = 2.727;
ro2 = 1744.92;
rp2 = 10.62;

gm3 = 0.1525;
ro3 = 23402.55;
rp3 = 572.9396;

% 10x10 matrix of system equations

% [vb1, ve1, vb2, ib1, ic1, ib2, ic2, ib3, ic3, id]

f_ac = [
    % F1
    [1/rp1, -1/rp1, 0, -1, 0, 0, 0, 0, 0, 0]; 
    % F2
    [gm1, -gm1 - 1/ro1, 0, 0, -1, 0, 0, 0, 0, 0];
    % f3
    [0, 1/rp2, -1/rp2, 0, 0, -1, 0, 0, 0, 0];
    % F4
    [0, gm2 + 1/ro2, -gm2, 0, 0, 0, -1, 0, 0, 0];
    % F5
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0];
    % F6
    [0, 0, -1/ro3, 0, 0, 0, 0, 0, 1, 0];
    % F7
    [1/1e3, -1/1e3, 0, 1, 0, 0, 0, 0, 0, 1];
    % F8
    [0, 0, 0, 0, 0, 1, 0, 0, -1, 1];
    % F9
    [1/1e3, -1/25 - 1/1e3, 0, 1, 1, -1, -1, 0, 0, 0];
    % F10
    [0, 1/100e3, 0, 0, 0, 0, 0, -1, 0, 0];
];


% Define solution vector to system f_AC
f_ac_sol = @(v_in) [
    0; % f1 O
    % -9/ro1; % f2
    0; % f2
    0; % f3
    0; % f4
    v_in/rp3; % f5
    gm3*v_in; % f6
    0; % f7
    0; % f8
    % -9/25; % f9
    0; % f9
    v_in/100e3; % f10
];

% Define input signal of 0.4mA amplitude and 10kHz frequency

v_in = 0.4e-3; % V
f_in = 10e3; % Hz
t = linspace(0, 3/f_in, 1000); % Time vector for three cycles

input = v_in * sin(2 * pi * f_in * t); % V
input = input(:); % Convert to column vector

% Gaussian Elimination for AC analysis
fprintf('--------------------------------------------------------------\n');
fprintf('Starting AC analysis using Gaussian elimination...\n');
fprintf('--------------------------------------------------------------\n');

% Initialize output voltage vector
x_ac_voltage = zeros(length(input), 1); 

% Initialize output current vector
x_ac_current = zeros(length(input), 1);

% Initialize input current vector
x_ac_input_current = zeros(length(input), 1);

% Loop through each value of the input signal and solve the system of equations using Gaussian elimination
for i = 1:length(input)
    x_ac_val = GaussianElimination(f_ac, f_ac_sol(input(i)));
    % x_ac_voltage(i) = 9 - x_ac_val(2); % Output voltage is 9 - VE1
    x_ac_voltage(i) = x_ac_val(2); % Output voltage is  VE1
    x_ac_input_current(i) = x_ac_val(8); % Input current is VB3
end

% Compute instantaneous output power across the 25-ohm resistor
input_power = (input.*x_ac_input_current);
output_power = (x_ac_voltage.^2) / 25;  % in Watts

% Plotting the Output Voltage
figure;
subplot(2,1,1);
hold on;
plot(t, x_ac_voltage);
plot(t, input);
xlabel('Time (s)');
ylabel('Output Voltage (V)');
title('Output Voltage vs. Time');

% Plotting the Output Power and input power vs time
subplot(2,1,2);
plot(t, output_power);
xlabel('Time (s)');
ylabel('Output Power (W)');
title('Output Power vs. Time');


figure
plot(t, input_power, 'r', 'DisplayName', 'Input Voltage');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Output Voltage vs. Time');
grid on;  

amp_in = max(input_power) - min(input_power);
amp_out = max(output_power) - min(output_power);

disp(size(input_power))
disp(size(output_power))

disp(amp_out/amp_in)

disp('--------------------------------------------------------------')
disp('Starting AC Solution using Gauss Seidel Iteration')
disp('--------------------------------------------------------------')

% Initialize output voltage vector
x_ac_voltage_GSI = zeros(length(input), 1); 

% Initialize output current vector
x_ac_current_GSI = zeros(length(input), 1);

% Initialize input current vector
x_ac_input_current_GSI = zeros(length(input), 1);

GSI_guess = x_ac_val; % Initial guess for Gauss-Seidel

% Loop through each value of the input signal and solve the system of equations using Gauss Seidel
for i = 1:length(input)
    [x_ac_val_GSI, iter_GSI] = GaussSeidel(f_ac, f_ac_sol(input(i)), 1e-9, 1e6, GSI_guess);
    fprintf('Iteration %d:\n', i);
    fprintf('Iterations for Gauss-Seidel = %d \n', iter_GSI);
    x_ac_voltage_GSI(i) = x_ac_val_GSI(2); % Output voltage is VE1
    disp('Output voltage from Gauss-Seidel:');
    disp(x_ac_voltage_GSI(i))
    disp('Output voltage from Gaussian elimination:');
    disp(x_ac_voltage(i))
    x_ac_input_current_GSI(i) = x_ac_val_GSI(8); % Input current is VB3
    GSI_guess = x_ac_val_GSI; % Update guess for next iteration
end

% Compute instantaneous output power across the 25-ohm resistor
input_power_GSI = (input.*x_ac_input_current_GSI);
output_power_GSI = (x_ac_voltage_GSI.^2) / 25;  % in Watts

% Compare values from GSI and Gaussian elimination for output voltage
figure
subplot(3,1,1);
plot(t, x_ac_voltage, 'r', 'DisplayName', 'Output Voltage (Gaussian Elimination)');
hold on;
plot(t, x_ac_voltage_GSI, 'b', 'DisplayName', 'Output Voltage (Gauss-Seidel)');
xlabel('Time (s)');
ylabel('Output Voltage (V)');
title('Output Voltage vs. Time');
grid on;
legend('Location', 'best');
subplot(3,1,2);
plot(t, input_power, 'r', 'DisplayName', 'Input Power (Gaussian Elimination)');
hold on;
plot(t, input_power_GSI, 'b', 'DisplayName', 'Input Power (Gauss-Seidel)');
xlabel('Time (s)');
ylabel('Input Power (W)');
title('Input Power vs. Time');
grid on;
legend('Location', 'best');
subplot(3,1,3);
plot(t, output_power, 'r', 'DisplayName', 'Output Power (Gaussian Elimination)');
hold on;
plot(t, output_power_GSI, 'b', 'DisplayName', 'Output Power (Gauss-Seidel)');
xlabel('Time (s)');
ylabel('Output Power (W)');
title('Output Power vs. Time');
grid on;
legend('Location', 'best');

% Graph output voltage from Gauss Seidel
figure
plot(t, x_ac_voltage_GSI, 'b', 'DisplayName', 'Output Voltage (Gauss-Seidel)');
xlabel('Time (s)');
ylabel('Output Voltage (V)');
title('Output Voltage vs. Time');
grid on;
legend('Location', 'best');

% Graph output voltage from Gaussian elimination
figure
plot(t, x_ac_voltage, 'r', 'DisplayName', 'Output Voltage (Gaussian Elimination)');
xlabel('Time (s)');
ylabel('Output Voltage (V)');
title('Output Voltage vs. Time');
grid on;
legend('Location', 'best');



amp_in_GSI = max(input_power_GSI) - min(input_power_GSI);
amp_out_GSI = max(output_power_GSI) - min(output_power_GSI);

fprintf('Solution from Gauss-Seidel:\n');
disp(amp_out_GSI/amp_in_GSI)