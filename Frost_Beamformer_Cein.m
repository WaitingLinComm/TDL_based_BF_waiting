%% Frost's beamformer
% from善晨
function [W,Y] = Frost_Beamformer_Cein(X,C,g,mu, w_opt)       % x = input signal, C = constraint matrix, g = constraint vector, mu = stepsize parameter

    L = numel(X(1,:));          % number of signal length
    w = C*inv(C'*C+0.0001*eye(numel(C(1,:))))*g;            % initialized weight vector
    MJ = numel(C(:,1));         % M*J (number of sensors * number of taps)
    M = numel(X(:,1));          % number of sensors
    J = MJ/M;                   % number of taps
    P = eye(MJ) - C*inv(C'*C+0.0001*eye(numel(C(1,:))))*C';   % alternative matrix
    Y = [];                     % output of beamformer
    W = [w];                    % weight vector for any time
    
    for i = 1:L-J+1
        x_n = reshape(flip(X(:,i:J+i-1),2),[MJ,1]);     % choose input signal iteratively
        y = w'*x_n;                             % output of beamformer
        Y = [Y y];                              % 
        e = x_n'*w;                             % error
        w = C*inv(C'*C+0.0001*eye(numel(C(1,:))))*g + P*(w-mu*e*x_n);       % update weight vector recursively     
        W = [W w];
        epsilon(i) = norm(w_opt - w);
    end
end