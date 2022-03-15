%% Covariance_MT_Calc
function [Rx, InvRx, InvRx_hat] = Correlation_MT_Calc_1(X_frames, M, J, Nframes, delta) 
% Assume source signal, interference, noise are correlated.
% The function calculate covariance matrix of the input data.
% 
% [Input]
% X_frames :input data with size: M * J * Nframes which is a three dimension array.
% M        :number of sensor.
% J        :number of taps.
% Nframes  :number of samples(Total observe time).  
% delta    :parameters of diagonal loading.
%
% [Output]
% Rx       :covariance matrix of the input data.
% InvRx    :inverse of covariance matrix of the input data 
% InvRx_hat:inverse of covariance matrix of the input data with diagonal loading

Rx_sum = zeros(M * J, M * J);
for i = 1 : Nframes            % Note: the TDL is filled with data after J + 1 time delay. 
    X = X_frames(:, :, i) ;    % Data matrix
    x = reshape(X, M * J, 1);  % Data vector.
    Rx_sum = Rx_sum +  x * x';
end
Rx = (1/Nframes) * Rx_sum;
InvRx = inv(Rx);
InvRx_hat = InvRx + delta * eye(M * J); % Diagonal loading

end


