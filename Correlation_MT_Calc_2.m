%% Covariance_MT_Calc
function [Rx, InvRx, InvRx_hat] = Correlation_MT_Calc_2(f_source, f_interference, E_s, E_v, sigma_eta, theta_s, theta_v, M, J, delta, Ts, mu) 
% Assume source signal, interference, noise are uncorrelated.
% The function calculate covariance matrix of the input data.
% 
% [Input]
% M        :number of sensor.
% J        :number of taps.
% delta    :parameters of diagonal loading.
% mu = d/cTs
%
% [Output]
% Rx       :covariance matrix of the input data.
% InvRx    :inverse of covariance matrix of the input data 
% InvRx_hat:inverse of covariance matrix of the input data with diagonal loading

R_s = zeros(M * J, M * J);
for i = 1:length(f_source)
    [~, a_s] = Steering_MT(M, J, f_source(i), theta_s, mu, Ts); % a_s is source steering vector(ULA)
    R_s = R_s + E_s^2 * (a_s * a_s');                           % Source signal autocorrelation matrix.
end
R_v = zeros(M * J, M * J);
for i = 1 : length(f_interference)
    [~, a_v] = Steering_MT(M, J, f_interference(i), theta_v(i), mu, Ts); % a_v1_1 is interference steering vector(ULA)
    R_v = R_v + E_v^2  * (a_v * a_v');                                   % Interference autocorrelation matrix.
end
R_eta = sigma_eta^2 * eye(M * J); % Noise autocorrelation matrix.
Rx = R_s + R_v + R_eta;
InvRx = inv(Rx);
InvRx_hat = InvRx + delta * eye(M * J); % with diagonal loading.

end


