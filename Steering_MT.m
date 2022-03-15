%% Covariance_MT_Calc
function [A, a] = Steering_MT(M, J, Omega, DOA, mu, Ts) 
% The function calculate steering matrix and steering vector of the TDL-based beamformer.
% 
% [Input]
% M        :number of sensor.
% J        :number of taps.
% DOA      :Direction of arrival.
% mu = d/cTs
%
% [Output]
% A        :steering matrix, size: (M, J) 
% a        :steering vector, size: (M * J, 1)

A = zeros(M, J);                    % Steering matrix
a = zeros(M * J, 1);                % Steering vector

tau0 = 0; 
% tau0 = -(J-1) * Ts/2; % If the middle point of the array is chosen as the zero-phase reference point. 
for k = 0:J-1
    A(:, k+1) = exp(-1i * Omega * ( tau0/Ts  + (0:M-1).' * mu * sin(DOA * pi / 180) + k));
end
a = reshape(A, M * J, 1);  % Concatenated steering vector

end