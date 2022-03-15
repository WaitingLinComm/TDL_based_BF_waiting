% From 善晨
clc ;
clear all;

%% Parameters
M = 10;         % number of sensor
J = 20;         % number of tap 
L = 200;       % data length

% desired signal
theta_s = 0;        % theta (degree)
freq_s_1 = 600;     % desired signal's frequency 1 (Hz)
freq_s_2 = 500;     % desired signal's frequency 2 (Hz)
% interference 1
theta_i_1 = 20;     % theta (degree)
freq_i_1_1 = 600;   % interference 1's frequency 1 (Hz)
freq_i_1_2 = 400;   % interference 1's frequency 2 (Hz)
% interference 2
theta_i_2 = -30;    % theta (degree)
freq_i_2_1 = 500;   % interference 2's frequency 1 (Hz)
freq_i_2_2 = 400;   % interference 2's frequency 2 (Hz)

freq_vec = [freq_s_1 freq_s_2 freq_i_1_1 freq_i_1_2 freq_i_2_1 freq_i_2_2];
freq_max = max(freq_vec);      % max frequency
T_s = 1/(2*freq_max);

%% received signal
s_1 = Data_Generation_single_tone(M,L,theta_s,freq_s_1,T_s);
s_2 = Data_Generation_single_tone(M,L,theta_s,freq_s_2,T_s);

i_1_1 = Data_Generation_single_tone(M,L,theta_i_1,freq_i_1_1,T_s);
i_1_2 = Data_Generation_single_tone(M,L,theta_i_1,freq_i_1_2,T_s);

% i_2_1 = Data_Generation_single_tone(M,L,theta_i_2,freq_i_2_1,T_s);
% i_2_2 = Data_Generation_single_tone(M,L,theta_i_2,freq_i_2_2,T_s);

% x = (sqrt(10)*s_1+sqrt(10)*s_2) + (10*i_1_1+10*i_1_2) + (10*i_2_1+10*i_2_2) + (1/sqrt(2)*randn(M,L+1)+1/sqrt(2)*randn(M,L+1));  %received signal
 x = (sqrt(10)*s_1+sqrt(10)*s_2) + (10*i_1_1+10*i_1_2) + (1/sqrt(2)*randn(M,L+1)+1/sqrt(2)*randn(M,L+1));  %received signal

save('test.mat','x');

%% Plot the Instantaneous Signal Power.
% This plots the instantaneous power for every element (M waveforms).
figure;
subplot(3, 1, 1);
plot((0:L), 10 * log10(abs(s_1).^ 2));  
grid on;
title('Instantaneous Source Signal Power');
xlabel('Sample Number');
ylabel('Output Power (dB)');

subplot(3, 1, 3);
plot((0:L), 10 * log10(abs(x).^ 2));  
grid on;
title('Instantaneous Source Signal Power');
xlabel('Sample Number');
ylabel('Output Power (dB)');
%% Generate single tone signal
function data_mat = Data_Generation_single_tone(M,L,theta,freq,T_s)        % M = number of sensor, L = signal length 
                                                                           % theta = DOA (degree), freq = frequency of single tone signal
                                                                           % T_s = sampling period
    DOA = abs(theta);   % make sure the DOA is positive. If theta is negative, we can reverse the result of the loop.
    data_mat = [];
    for m = 0:M-1
        tau = m*sin(DOA*pi/180)*T_s;       % spatial delay of mth antenna 
        n = ceil(m*sin(DOA*pi/180));       % number of zero padding for no received signal
        received_L = linspace(n*T_s-tau, (n*T_s-tau)+(L-n)*T_s, L-n+1);  % sampling length when we received the signal
        signal = [zeros(1,n)  exp(j*2*pi*freq*received_L)];   % signal received by mth antenna
        data_mat = [data_mat; signal];
    end
        
    if theta<0
        data_mat = flip(data_mat,1);
    end
    
end


