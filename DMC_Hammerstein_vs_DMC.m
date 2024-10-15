%% Nonlinear Hammerstein Dynamic Matrix Control (NLH-DMC)
% Author: Sergio Andres Castaño Giraldo
% YouTube: https://www.youtube.com/@SergioACGiraldo
% Website: https://controlautomaticoeducacion.com/
% GitHub: https://github.com/sergioacg
%
% This MATLAB script implements a Nonlinear Hammerstein Dynamic Matrix Control (NLH-DMC)
% for a tank level control system. The controller uses a nonlinear static gain model (Hammerstein model)
% combined with a linear dynamic model in a predictive control setup.
% It simulates the control of tank level by adjusting the valve opening based on predictive control theory.

clc;
clear all;
close all;

% Define tank system parameters
k1 = 0.05;    % Flow coefficient for valve 1
k2 = 0.015;   % Flow coefficient for valve 2
g = 10;       % Gravity constant
A = 0.5;      % Cross-sectional area of the tank
alpha = 10;   % Proportionality constant for the valve

% Pack parameters into a structure for easier passing to functions
params.k1 = k1;
params.k2 = k2;
params.g = g;
params.A = A;
params.alpha = alpha;

% Sampling time
Ts = 10;
% Polynomial coefficients for the nonlinear valve gain (Hammerstein model)
polynomial = [4.7269, -6.4556, 4.6015, -0.8817, 0.2077, 0.0208];
%polynomial = [-16.8154   32.6937  -19.8106    5.1759   -0.2969    0.0259];

% Linear dynamic model of the system as a transfer function
Gp = tf(1, [17.69, 1]); %Hammerstein
Gdmc = tf(1.769, [17.69, 1]); %Classical DMC

% Convert continuous model to discrete with zero-order hold
ftz = c2d(Gp, Ts);
ftz_dmc = c2d(Gdmc, Ts);

% Delay inherent in the model
d = ftz.iodelay;

%% Predictive Control Setup

% Prediction horizon (typically a few times the system's time constant)
P = 10;                   
% Control horizon (often less than the prediction horizon)
N = 5;                    
% Error weighting matrix
delta = 1 * eye(P);         
% Control effort weighting factor
lamda = 0.8;                

Qd = eye(P) * delta;
Ql = eye(N) * lamda;

% Number of samples to use from the step response
Nm = P + d;  

%% Gain Control
[K1, gi] = GainControlDMC(ftz, P, N, Qd, Ql)


%% Control Loop

% Comparative Analysis of Classical DMC and NLH-DMC
for j = 1: 2
    
    if j == 2
        P = 10;
        N = 5;
        delta = 1 * eye(P);
        lamda = 2;   
        Qd = eye(P) * delta;
        Ql = eye(N) * lamda;
        % Gain Control
        [K1, gi] = GainControlDMC(ftz_dmc, P, N, Qd, Ql);
    end

    % Initialization of variables
    nit = 100;  % Number of iterations
    operationPoint = 0.6;  % Initial valve opening
    a2 = 0.5;
    u = [operationPoint a2];
    t = 0:Ts:10;
    equilibriumHeight = fsolve(@(h) tank_model(t, h, u, params), 0.01);

    u(1:nit) = operationPoint; % Initialize control input

    initialTankHeight = equilibriumHeight; % Initial tank height
    currentHeight = initialTankHeight;
    ym = currentHeight * ones(1, nit); % Tank height measurements

    % Reference signal setup
    r(1:nit) = initialTankHeight;  % Set initial reference
    r(6:nit) = initialTankHeight + 0.4; 
    r(50:nit) = initialTankHeight - 0.1;

    xt(1:nit) = initialTankHeight; % Initialize pseudo output

    dx = zeros(1, Nm); % Control action (Delta U Free)

    a2(1:nit) = 0.5;     % Constant inlet valve opening (considered as a disturbance)
    a2(35:nit) = 0.55;
    a2(80:nit) = 0.4;   % Change in disturbance

    inc_x(1:nit) = 0;

    w=0;  %For future references
    for k = 2:nit
        % Simulate the Nonlinear Model
        uin = [u(k-1), a2(k)];
        tspan = [0 Ts]; % Sampling period to integrate
        [ts, x] = ode45(@(t, x) tank_model(t, x, uin, params), tspan, currentHeight);
        % Update the current height and apply saturation
        currentHeight = max(0, min(x(end, :), 1));
        %currentHeight = x(end, :);
        ym(k) = currentHeight; % Update measured tank height

        
         %% Calculate Free Response
         f = zeros(P, 1); % Initialize free response vector

         for kk = 1:P
            vect_g = zeros(Nm, 1);  % Initialize vect_g for each kk

            % Calculate the differences within the range of gi
            for i = 1:Nm-P
                vect_g(i) = gi(kk+i) - gi(i);
            end
            % Handle the remaining elements where kk+i exceeds the length of gi
            for i = Nm-P+1:Nm
                vect_g(i) = gi(Nm) - gi(i);
            end

            % Calculate the free response by using the dot product of vect_g and dx
            f(kk) = ym(k) + vect_g' * dx';  
         end


        % Calculate Control Increment
        % No future reference available, thus using only current reference
        inc_x(k) = K1 * (r(k + w) * ones(P, 1) - f);  % Corrected to use column vector for ones

        % Update the control action history
        aux_2 = dx(1:Nm-1);  % Capture all but the last element of dx
        dx = [inc_x(k) aux_2];  % Prepend the new control increment to the history

        % Compute the Pseudo-Output xt
        if k == 1
            xt(k) = inc_x(k);  % Initialize xt at the first iteration
        else
            xt(k) = xt(k-1) + inc_x(k);  % Accumulate control action to pseudo-output
        end

         % Calculate and apply control action
         xp = polynomial;
         xp(end) = polynomial(end) - xt(k);  % Adjust polynomial for current control

         % Calculate the roots of the polynomial
        ur = roots(xp);  

        % Find all real roots
        real_roots = ur(imag(ur) == 0);  % Only real roots

        % Find real roots within the range [0, 1]
        ww = find(real_roots >= 0 & real_roots <= 1);  

        if ~isempty(ww)  % Check if there is at least one suitable root within the range
            % If multiple roots within range, you might choose the smallest, largest, etc.
            selected_root = min(real_roots(ww));  % Example: Choose the smallest
        else
            % If no roots are within the range, find the real root closest to the range
            if isempty(real_roots)
                selected_root = NaN;  % No real roots at all
            else
                % Find the closest real root to the interval [0, 1]
                [~, idx] = min(abs(real_roots - 0.5));  % Closest to the middle of the range
                selected_root = real_roots(idx);
            end
        end

         
         
         if j == 1
            u(k) = selected_root;  % Apply the real root as control action
         else
             u(k) =  u(k-1) + inc_x(k);
         end
         % Control saturation
         if u(k) > 1
             u(k) = 1;
         elseif u(k) < 0
             u(k) = 0;
         end
    end
    if j == 1
        ym_ham = ym;
        u_ham = u;
        inc_xH = inc_x;
    else
        ym_classical = ym;
        u_classical = u;
    end
end

%% Plotting results for both control strategies
t = 0:Ts:(nit-1)*Ts;
figure;
subplot(2,1,1);
plot(t, r, '--k', t, ym_ham, '-r', t, ym_classical, '--b', 'Linewidth', 2);
xlabel('Time (s)');
ylabel('Tank Height');
legend('Setpoint', 'NLH-DMC', 'Classical DMC', 'Location', 'Northeast');
grid on;

subplot(2,1,2);
stairs(t, u_ham, 'r', 'Linewidth', 3);
hold on;  % Keep the plot to add the next control action series
stairs(t, u_classical, '--b', 'Linewidth', 3);
hold off;
xlabel('Time (s)');
ylabel('Control Input');
legend('Valve Position NLH-DMC', 'Valve Position Classical DMC');
grid on;

% figure
% plot(t, inc_xH, '-k', 'Linewidth', 2);
% xlabel('Time (s)');
% ylabel('Pseudo-Output Increment');
