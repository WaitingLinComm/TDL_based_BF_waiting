%% Workspace Initialization.
clc; clear; close all;
%% Array Parameters.
M = 10;    % Number of array sensors.
J = 20;    % Number of taps.
%% Signal Parameters.
Nsamples = 2000; % Number of samples(Total observe time).  

% Source signal
%f_source = [600 500]; % Hz. Frequency vector.
f_source = [600]; 
theta_s = 0;   % deg

% Interference
%f_interference = [600 400 300]; % Frequency vector.
%theta_v = 20 * ones(1, length(f_interference));
f_interference = [600];
theta_v = [20]; % deg. Note: length(theta_v) should be equal to length(f_interference).

% System parameters.
fmax = max([f_source f_interference]);       % Sampling frequency.
Ts   = 1/(2 * fmax);                         % Sampling period.

% Signal Energy.
SNR       = 10;        % Signal to Noise Ratio 
INR1      = 20;        % Interference to Noise Ratio 
E_s       = 10^(SNR/20);  % Energy/Voltage/Amplitude of source signal is E_s (= sigma_s). Power of source signal is P_s = Es^2 = sigma_s^2. 
E_v       = 10^(INR1/20); % Energy/Voltage/Amplitude of interference  is E_v (= sigma_v). Power of interference is P_v = E_v^2 = sigma_v^2.
sigma_eta = 1;      % E_eta = sigma_eta = 1; % P_eta = 1;

% Frequency information
Omega_L     = 0.25 * pi;   
Omega_H     = pi; 
Freq_range  = Omega_H - Omega_L;
N_Omega     = 100;                                          % Number of frequency bins to observe.
Omega_I     = Omega_L : Freq_range/(N_Omega - 1) : Omega_H; % The frequency range of interest. 
Omega_r     = 0.25 * pi;                                    % Reference frequency.
mu          = 1;                                            % mu = d/cTs
%% Generate wideband signal.
% Generate wideband source signal.
data_S = zeros(M, Nsamples);
for i = 1 : length(f_source)
    data_S = data_S + E_s * Single_tone_data_generator(M, Nsamples, theta_s, f_source(i), Ts); 
end
% Generate wideband interference.
data_V = zeros(M, Nsamples);
for i = 1 : length(f_interference)
    data_V = data_V + E_v * Single_tone_data_generator(M, Nsamples, theta_v(i), f_interference(i), Ts);
end
% Noise: uncorrelated unit power thermal noise samples with a complex Gaussian distribution.
data_Eta = (sigma_eta/sqrt(2)) * (randn(M, Nsamples) + 1i * randn(M, Nsamples)); 

data_X = data_S + data_V +  data_Eta;
% Save data into .mat file/ load data from .mat file.
%{
% save('test_Data.mat','data_X');
% Load data(data stream)
% data = load('test_Data.mat'); 
% data_X = data.data_X; 
% if (M ~= size(data_X, 1))
%     printf("Warning: Loaded data is not correspond to the number of array sensors!");
% end
%}
% Generate data frame 
Nsamples    = size(data_X, 2);                      % Nsamples = 100; % Number of samples(Total observe time).  
Nframes     = Nsamples - 3 * J + 3;                 % Number of frames.
[X_frames]  = Data_Frame_GEN(data_X, Nframes, M, J);% size: M * J * Nframes which is a three dimension array.

[S_frames]  = Data_Frame_GEN(data_S, Nframes, M, J);
[V_frames]  = Data_Frame_GEN(data_V, Nframes, M, J);
[Eta_frames]= Data_Frame_GEN(data_Eta, Nframes, M, J);
%% Calculate correlation matrix.
% Rx        :covariance matrix of the input data.
% InvRx     :inverse of covariance matrix of the input data 
% InvRx_hat :inverse of covariance matrix of the input data with diagonal loading
delta = 0.00001; % parameters of diagonal loading.

% Sol1: Assume source signal, interference, noise are not uncorrelated.
[Rx, InvRx, InvRx_hat] = Correlation_MT_Calc_1(X_frames, M, J, Nframes, delta); 
[Rs, ~, ~]   = Correlation_MT_Calc_1(S_frames, M, J, Nframes, delta); 
[Rv, ~, ~]   = Correlation_MT_Calc_1(V_frames, M, J, Nframes, delta); 
[Reta, ~, ~] = Correlation_MT_Calc_1(Eta_frames, M, J, Nframes, delta); 
%Reta = sigma_eta ^2 * eye(M * J);

% So12: Assume source signal, interference, noise are uncorrelated.
%[Rx, InvRx, InvRx_hat] = Correlation_MT_Calc_2(f_source, f_interference, E_s, E_v, sigma_eta, theta_s, theta_v, M, J, delta, Ts, mu); 
%% Frost Beamformer(Frost Beamformer is only suitable for wideband input signal with DOA = 0.)
w_opt = LCMV_Beamformer(InvRx_hat, M, J); 
%% Freq domain Beamformer
% Only contains source signal in constraint matrix
% Constraint matrix
C = zeros(M * J, N_Omega);
for r = 1 : N_Omega
    [A_s, a_s] = Steering_MT(M, J, Omega_I(r), theta_s, mu, Ts); 
    C(:, r)    = a_s;
end
nu         = 0.000004;
iterations = Nframes;           % Number of input data frame.
g          = ones(N_Omega, 1);  % size:(N_Omega, 1) 
[w_frost, outputpower, epsilon, SINR_buffer]  = Frost_Beamformer(w_opt, C, g, M ,J ,nu, iterations, N_Omega, delta, Rx, Rs, Rv, Reta, X_frames, S_frames, V_frames, Eta_frames); 

figure;
plot(outputpower, 'r.');
figure;
plot(epsilon, 'b.'); 
figure;
plot(10 * log10(SINR_buffer), 'b.'); 

% from善晨
%{
g       = ones(N_Omega, 1);         % size:(N_Omega, 1) 
nu      = 0.000004;
[W_test, Y_test] = Frost_Beamformer_Cein(data_X, C, g, nu, w_opt); 
w_opt_test = W_test(:,end);
%}

% Consider source signal +  one interference in constraint matrix.
%{
C = zeros(M * J, N_Omega * 2);
for r = 1 : N_Omega
    [A_s, a_s] = Steering_MT(M, J, Omega_I(r), theta_s, mu, Ts); 
    C(:, r) = a_s;
end
for r = 1 : N_Omega
    [A_s, a_v] = Steering_MT(M, J, Omega_I(r), theta_v, mu, Ts); 
    C(:, N_Omega + r) = a_v;
end
nu         = 0.000004;
iterations = Nframes;                                % Number of input data frame.
g          = [ones(N_Omega, 1); zeros(N_Omega, 1)];  % size:(2* N_Omega, 1)
[w_frost, outputpower, epsilon, SINR_buffer]  = Frost_Beamformer2(w_opt, C, g, M ,J ,nu, iterations, N_Omega, delta, Rx, Rs, Rv, Reta, X_frames, S_frames, V_frames, Eta_frames); 
figure;
plot(outputpower, 'r.');
figure;
plot(epsilon, 'b.'); 
figure;
plot(10 * log10(SINR_buffer), 'b.'); 
%}

%% Response Variation(RV) Beamformer
% beta = 10; % 1. 0.1
% w_RV = RV_Beamformer(beta, Omega_I, Omega_r, theta_s, Ts, mu, Rx, M, J); 
%% Beampattern(Frost/RV Beamformer)
theta_grid = (-90 : 0.1 : 90);      % Plotted theta_grid interval(in degrees).
Ntheta = length(theta_grid);        % Ntheta = 1801
A = zeros(M, J);                    % Steering matrix
a = zeros(M * J, 1);                % Steering vector
Beam_Pat_LCMV = zeros(N_Omega, Ntheta);
Beam_Pat_frost  = zeros(N_Omega, Ntheta);
Beam_Pat_RV  = zeros(N_Omega, Ntheta);
for i = 1:Ntheta
    for r = 1:N_Omega
        [A, a] = Steering_MT(M, J, Omega_I(r), theta_grid(i), mu, Ts); 
        Beam_Pat_LCMV(r, i) = w_opt' * a;
        Beam_Pat_frost(r, i) = w_frost' * a;
        %Beam_Pat_RV(r, i) = w_RV' * a;
    end
end
%% Plot Beampattern.
%%%%%%%%%%%%%%%% Frost beamformer %%%%%%%%%%%%%%%%

figure
mesh(theta_grid, Omega_I, 10 * log10(abs(Beam_Pat_frost).^2));
title('Beampattern of Frost beamformer');
ylabel('Normalized frequency: $\Omega/\pi$','interpreter','latex');
xlabel('$DOA: \theta$ (deg)','interpreter','latex');
zlabel('$|B(\Omega, \theta)|(dB)$','interpreter','latex');

xlim([-90 90]);
set(gca,'xtick',[-90 -45 0 45 90]); 
% set(gca,'xticklabel',{'-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2'}); 

ylim([Omega_L Omega_H]);
set(gca,'ytick',[Omega_L, (Omega_L + Omega_H)/2, Omega_H]); 
set(gca,'yticklabel',{'0.25', '0.5' ,'1'}); 

%%%%%%%%%%%%%%%% LCMV beamformer %%%%%%%%%%%%%%%%
%{
figure
mesh(theta_grid, Omega_I, 10 * log10(abs(Beam_Pat_LCMV).^2));
title('Beampattern of LCMV beamformer');
ylabel('Normalized frequency: $\Omega/\pi$','interpreter','latex');
xlabel('$DOA: \theta$ (deg)','interpreter','latex');
zlabel('$|B(\Omega, \theta)|(dB)$','interpreter','latex');

xlim([-90 90]);
set(gca,'xtick',[-90 -45 0 45 90]); 
% set(gca,'xticklabel',{'-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2'}); 

ylim([Omega_L Omega_H]);
set(gca,'ytick',[Omega_L, (Omega_L + Omega_H)/2, Omega_H]); 
set(gca,'yticklabel',{'0.25', '0.5' ,'1'}); 
%}
%%%%%%%%%%%%%%%% RV beamformer %%%%%%%%%%%%%%%%
%{
figure
mesh(theta_grid, Omega_I, 10 * log10(abs(Beam_Pat_RV).^2));
title('Beampattern of RV beamformer');
ylabel('Normalized frequency: $\Omega/\pi$','interpreter','latex');
xlabel('$DOA: \theta$ (deg)','interpreter','latex');
zlabel('$|B(\Omega, \theta)|(dB)$','interpreter','latex');

xlim([-90 90]);
set(gca,'xtick',[-90 -45 0 45 90]); 
% set(gca,'xticklabel',{'-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2'}); 

ylim([Omega_L Omega_H]);
set(gca,'ytick',[Omega_L, (Omega_L + Omega_H)/2, Omega_H]); 
set(gca,'yticklabel',{'0.25', '0.5' ,'1'}); 
%}
