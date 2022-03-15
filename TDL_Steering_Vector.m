%% Generate single tone signal
function [a] = TDL_Steering_Vector(M, J, f, fmax)        % M = number of sensor, L = signal length                                                              % theta = DOA (degree), freq = frequency of single tone signal
                                                                           % T_s = sampling period
    A = zeros(M, J);
    Omega_s = pi * (f/fmax);
    for k = 0:J-1
        A(:, k+1) = exp(-1i * Omega_s * tau0/Ts) * exp(-1i * Omega_s * (0:M-1).' * mu * sin(theta_s)) * exp(-1i * Omega_s * k);
    end
    a = reshape(A, M * J, 1);
end
