function dh = tank_model(t, h, u, par)
    % Load parameters
    k1 = par.k1; % Flow coefficients for the valves
    k2 = par.k2; % Flow coefficients for the valves
    A = par.A;   % Cross-sectional area of the tank
    g = par.g;   % Gravitational acceleration
    alpha = par.alpha; % Proportionality parameter for the equal percentage valve

    % Model inputs
    a1 = alpha ^ (u(1) - 1); % Equal percentage opening of the valve
    a2 = u(2);               % Linear opening of the inlet valve

    % Model equation
    dh = (1 / A) * (a1 * k1 - a2 * k2 * sqrt(2 * g * h));
end
