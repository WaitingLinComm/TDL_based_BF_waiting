function [w_frost, outputpower, epsilon, SINR_buffer] = Frost_Beamformer2(w_opt, C, g, M ,J ,nu, iterations, N_Omega, delta, Rx, Rs, Rv, Reta, X_frames, S_frames, V_frames, Eta_frames) 
% The function implement LMS-type solution of Frost beamformer.
% 
% [Input]
% M        :number of sensor.
% J        :number of taps. 
% delta    :parameters of diagonal loading.
%
% [Output]
% w_frost    :output weight.


W       = zeros(M * J, iterations); % size:(M, iterations) 
w_0     = C * inv(C' * C + delta * eye(N_Omega * 2)) * g; % size:(M * J, 1).  % with diagonal loading. 
% w_0     = C * inv(C' * C ) * g;   % size:(M * J, 1) % without diagonal loading.
W(:,1)  = w_0;

P       = eye(M * J) - C * inv(C' * C + delta * eye(N_Omega * 2)) * C'; % size:(M * J, M * J) 
% P       = eye(M * J) - C * inv(C' * C) * C'; % size:(M * J, M * J) 

e           = zeros(1, iterations); 
y_buffer    = zeros(1, iterations); 
epsilon     = zeros(1, iterations); 
SINR_buffer = zeros(1, iterations);
outputpower = zeros(1, iterations); s_power   = zeros(1, iterations);
v_power     = zeros(1, iterations); eta_power = zeros(1, iterations);

for i = 1 : iterations
    X = X_frames(:, :, i); x = reshape(X, M * J, 1);
    y_buffer(i) = W(:, i)' * x;
    e(i) = x' * W(:, i);
    W(:, i + 1) = w_0 + P * (W(:, i) - nu * e(i) * x);
    
    % outputpower(i) = norm( W(:, i)' * Rx * W(:, i)); % outputpower(i) =  (W(:, i)' * x) * conj(W(:, i)' * x);
    outputpower(i) = norm( W(:, i)' * x);  
    epsilon(i) = norm(w_opt - W(:, i));
    
    % Calculate SINR
    S = S_frames(:, :, i); s = reshape(S, M * J, 1);
    V = V_frames(:, :, i); v = reshape(V, M * J, 1);
    Eta = Eta_frames(:, :, i); eta = reshape(Eta, M * J, 1);

    % Sol1:
    s_power(i) = norm( W(:, i)' * s);  
    v_power(i) = norm( W(:, i)' * v);  
    eta_power(i) = norm( W(:, i)' * eta); 
    SINR_buffer(i) = s_power(i) / (v_power(i) + eta_power(i));
    % Sol12:
    %{
    s_power(i) = norm( W(:, i)' * Rs * W(:, i));  
    v_power(i) = norm( W(:, i)' * Rv * W(:, i));  
    eta_power(i) = norm( W(:, i)' * Reta * W(:, i)); 
    SINR_buffer(i) = s_power(i) / (v_power(i) + eta_power(i));
    %}
    % Sol3:
    %{
    Rs = s * s';
    Rv = v * v';
    Reta = eta * eta';
    SINR_buffer(i) = W(:, i)' * Rs * W(:, i)/ (W(:, i)' * (Rv +  Reta)  * W(:, i));
    %}
end

w_frost = W(:, end);
% figure;
% plot(outputpower, 'r.');
% figure;
% plot(epsilon, 'b.'); 
% figure;
% plot(10 * log10(SINR_buffer), 'b.'); 

end


