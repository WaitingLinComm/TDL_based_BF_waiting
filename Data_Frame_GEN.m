%% Covariance_MT_Calc
function [X_frames] = Data_Frame_GEN(data_X, Nframes, M, J) 
% The function translates input data stream to data frame with size: M * J * Nframes.
% 
% [Input]
% data_X   :input data stream.
% Nframes  :number of frames.
% M        :number of sensor.
% J        :number of taps. 
%
% [Output]
% X_frames :data frame with size: M * J * Nframes which is a three dimension array.

X_frames = zeros(M, J, Nframes);

for i = 1 : Nframes % Note: the TDL is filled with data after J + 1 time delay.
    X_tmp = data_X(:, J + i : J + i + J - 1); % size: M * J
    X = fliplr(X_tmp); % Data matrix, size: M * J. % Note: the data matrix should be reversed.
    X_frames(:, :, i) = X;  % size: M * J * Nframes which is a three dimension array.
end

end


