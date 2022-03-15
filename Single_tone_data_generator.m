function data_matrix = Single_tone_data_generator(M, Nsamples, DOA, f, Ts)                                                                                   
% M: number of sensor.
% Nsamples: number of samples(Total observe time).  
% theta: DOA (degree).
% f: frequency of single tone signal.
% Ts = 1/fmax: sampling period.

data_matrix = zeros(M, Nsamples);                                             
signal_0 = exp(1i * 2 * pi * f * 1); % Single tone signal received by the first sensors at time index 0.
for n = 1 : Nsamples
    tau0 = sin(DOA * pi / 180) * Ts;
    a = exp(-1i * 2 * pi * f * (0:M-1).' * tau0);
    data_matrix(:,n) = a * exp(1i * 2 * pi * f * (n - 1) * Ts) * signal_0;
end   

end