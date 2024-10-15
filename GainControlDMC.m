function [K1, gi] = GainControlDMC(ftz, P, N, Qd, Ql)
    
    Ts = ftz.Ts;
    d = ftz.iodelay;
    % Step response vector up to P+d time steps
    gi = step(ftz, (P + d)*Ts);  
          

    % Construct matrix G using the step response data
    G = zeros(P, N);
    G(:, 1) = gi(1 + d : P + d);  % Fill the first column with step response values

    for i = 2:N
        for j = 2:P
            G(j, i) = G(j - 1, i - 1);  % Propagate values diagonally to fill the matrix
        end
    end

    %% Calculate the gain matrix Mn
    % Unconstrained cost function optimization
    Mn = inv(G' * Qd * G + Ql) * G' * Qd; 

    % Control gain vector (first row of Mn)
    K1 = Mn(1, :);