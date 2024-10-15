%% Tank Level Control System Simulation
%
% Author: Sergio Andres Castaño Giraldo
% YouTube: https://www.youtube.com/@SergioACGiraldo
% Website: https://controlautomaticoeducacion.com/
% GitHub: https://github.com/sergioacg
%
% This MATLAB script simulates the dynamic behavior of a tank level control
% system using a nonlinear and linearized model. 
% The script calculates tank height based on different valve openings using 
% both a Hammerstein model and a traditional linearized approach.
% It first establishes a static relationship between the valve opening and 
% the tank height using a nonlinear model solved by fsolve. 
% Then, it linearizes the model at a specific operating point and compares 
% the dynamic response of the nonlinear model, 
% Hammerstein model, and linearized model under various test conditions to 
% illustrate the effects of nonlinearity and control strategies.


clc;
clear all;
close all;

% Tank variables
k1 = 0.05;    % Flow coefficient for valve 1
k2 = 0.015;   % Flow coefficient for valve 2
g = 10;       % Gravity constant
A = 0.5;      % Cross-sectional area of the tank
a1 = 0:0.05:1; % Opening of the inlet valve range
a2 = 0.5;     % Opening of the outlet valve
alpha = 10;   % Proportionality constant for the valve

% Structure to hold parameters
params.k1 = k1;
params.k2 = k2;
params.g = g;
params.A = A;
params.alpha = alpha;

% Static graph varying a1 vs tank height using fsolve to solve the model
initialHeight = 0; % initial height
heights = zeros(1, length(a1));
Ts = 0.5;
t = 0:Ts:10;
options = optimoptions('fsolve', 'Display', 'off');

for i = 1:length(a1)
    valveSettings = [a1(i) a2];
    height = fsolve(@(h) tank_model(t, h, valveSettings, params), initialHeight, options);
    heights(i) = max(0, min(height, 1));
end

% Find a fifth-order polynomial that fits the height data
polynomialCoefficients = polyfit(a1, heights, 5);
Xp = polyval(polynomialCoefficients, a1);

% Print the fitting polynomial
fprintf('Fitting polynomial: %f*a1^5 + %f*a1^4 + %f*a1^3 + %f*a1^2 + %f*a1 + %f\n', polynomialCoefficients);

figure(1)
plot(a1, heights, 'linewidth', 2)
hold on
plot(a1, Xp, '--r', 'linewidth', 2)
legend('Model', 'Polynomial fit')
xlabel('inlet valve opening')
ylabel('Tank height')
title('Input a1 vs. Tank height')

%% Linearize the model using Jacobian
operationPoint = 0.6; % Define an operation point
% Calculate an equilibrium height for the operation point
u = [operationPoint a2];
equilibriumHeight = fsolve(@(h) tank_model(t, h, u, params), initialHeight, options);

% Symbolic variables for tank height (x) and valves opening (u)
syms x u

% Function for Jacobian calculation
f1 = (1/params.A) * (u*params.k1 - a2*params.k2*sqrt(2*params.g*x));

% Calculation of the Jacobian
A = jacobian(f1, x);
B = jacobian(f1, u);

% Substitute the equilibrium points in the Jacobian
A = double(subs(A, [x, u], [equilibriumHeight, operationPoint]));
B = double(subs(B, [x, u], [equilibriumHeight, operationPoint]));

% Determine the Transfer Function
[num, den] = ss2tf(A, B, 1, 0);
G = tf(num/den(2), den/den(2));
Gham = G;
Gham.num = 1; % Set numerator to unity as only the dynamic part is of interest


% Simulate the nonlinear model (ode45) and compare with the linearized model at various points
nit = 1000; % Number of iterations
t = 0:Ts:(nit-1)*Ts; % Time vector
% Valve openings
a1(1:50) = operationPoint; 
a1(51:nit) = operationPoint + 0.02; 
a1(301:nit) = operationPoint + 0.2; 
a1(551:nit) = operationPoint - 0.4;
ulin = zeros(1, nit);

initialTankHeight = equilibriumHeight; % Initial tank height
currentHeight = initialTankHeight;
modelHeights = currentHeight * ones(1, nit); % Tank height from the model
heightsHam = currentHeight * ones(1, nit); % Tank height from the Hammerstein model
heightsLin = currentHeight * ones(1, nit); % Tank height from the linearized model

for k = 2:nit
    % Simulation of the Nonlinear Model
    u = [a1(k), a2];
    tspan = [0 Ts]; % Sampling period to integrate
    [ts, x] = ode45(@(t, x) tank_model(t, x, u, params), tspan, currentHeight);
    currentHeight = x(end, :); % Update the initial condition
    currentHeight = max(0, min(currentHeight, 1));
    modelHeights(k) = currentHeight; % The process output is the second state variable

    % Simulation of the Linear Model
    % Nonlinear gain of Hammerstein with the fitting polynomial
    if (a1(k) - operationPoint) ~= 0
        Ku = (polyval(polynomialCoefficients, a1(k)) - initialTankHeight) / (a1(k) - operationPoint);
    else
        Ku = 1;
    end
    % Calculate the difference in valve opening
    if (a1(k) - operationPoint) ~= 0
        ulin(k) = (a1(k) - operationPoint);
    end
    y = lsim(Ku * Gham, ulin(1:k), t(1:k), 'zoh') + initialTankHeight;
    heightsHam(k) = y(end);
    
    y = lsim(G, ulin(1:k), t(1:k), 'zoh') + initialTankHeight;
    heightsLin(k) = y(end);
end

figure(2)
plot(t, modelHeights, 'linewidth', 2)
hold on
plot(t, heightsHam, '--r', t, heightsLin, ':k', 'linewidth', 2)
legend('Nonlinear model', 'Hammerstein model', 'Linearized model')
xlabel('Time [s]')
ylabel('Tank height')
title('Comparison of nonlinear and linearized models')
