function [w_opt] = LCMV_Beamformer(InvRx, M, J) 
% The function calculate the optimal weight of Frost beamformer.
% Note: Frost Beamformer is only suitable for wideband input signal with DOA = 0.
% 
% [Input]
% InvRx    :inverse of covariance matrix of the input data with or without diagonal loading
% M        :number of sensor.
% J        :number of taps. 
% delta    :parameters of diagonal loading.
%
% [Output]
% w_opt    :optimal weight.

C = zeros(M * J, J);
for j = 1:J
    col = zeros(M * J, 1);
    col((j - 1) * M + 1: (j - 1) * M + M) = 1;
    C(:, j) = col;
end
f = zeros(J, 1);
f(1,1) = 1;
w_opt = InvRx * C * inv(C' * InvRx * C) * f;

end


