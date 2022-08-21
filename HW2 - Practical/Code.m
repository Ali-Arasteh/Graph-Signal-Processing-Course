gsp_start;
%% a)
%% 1)
A = [0 1 1 0 0 0 0 1; 1 0 1 1 1 0 0 1; 1 1 0 1 0 0 0 0; 0 1 1 0 1 1 0 1;
    0 1 0 1 0 1 1 1; 0 0 0 1 1 0 1 0; 0 0 0 0 1 1 0 0; 1 1 0 1 1 0 0 0];
W = A;
C = [0 0; 1 1; 0 2; 2 2; 3 1; 4 2; 4 1; 2 0];
G = gsp_graph(W, C);
G.plotting.edge_color = 'b';
G.plotting.vertex_color = 'k';
figure();
gsp_plot_graph(G);
title('G');
%% 2)
G = gsp_compute_fourier_basis(G);
x = 2 * G.U(:, 1) + G.U(:, 2);
%% 3)
SNR = 10;
signal_P = 1 / length(x) * (x') * x;
noise_variance = signal_P / 10^(SNR / 10);
noisy_x = x + randn(length(x), 1) * sqrt(noise_variance);
figure();
gsp_plot_signal(G, x);
title('x');
figure();
gsp_plot_signal(G, noisy_x);
title('noisy x');
figure();
stem(1:length(x), x);
ylim([-2 2])
title('x');
figure();
stem(1:length(noisy_x), noisy_x);
ylim([-2 2])
title('noisy x');
%% 4)
[VW,DW] = eig(W);
Wnorm = 1 / max(abs(DW), [], 'all') * W;
DWnorm = DW / max(abs(DW), [], 'all');
x_Wnorm_frequency = VW \ x;
figure();
stem(diag(DWnorm), x_Wnorm_frequency)
xlim([-1 1])
ylim([-4 4])
title('frequency spectrum of x (Wnorm)');
noisy_x_Wnorm_frequency = VW \ noisy_x;
figure();
stem(diag(DWnorm), noisy_x_Wnorm_frequency)
xlim([-1 1])
ylim([-4 4])
title('frequency spectrum of noisy x (Wnorm)');
x_L_frequency = G.U \ x;
figure();
stem(G.e, x_L_frequency)
xlim([0 8])
ylim([-2 2])
title('frequency spectrum of x (L)');
noisy_x_L_frequency = G.U \ noisy_x;
figure();
stem(G.e, noisy_x_L_frequency)
xlim([0 8])
ylim([-2 2])
title('frequency spectrum of noisy x (L)');
%% 5)
h_Wnorm_frequency = [0 0 0 0 0 0 1 1]';
figure();
stem(diag(DWnorm), h_Wnorm_frequency)
xlim([-1 1])
ylim([0 1])
title('frequency response of h (Wnorm)');
h_L_frequency = [1 1 0 0 0 0 0 0]';
figure();
stem(G.e, h_L_frequency)
xlim([0 8])
ylim([0 1])
title('frequency response of h (L)');
%% 6)
filtered_noisy_x_Wnorm_frequency = noisy_x_Wnorm_frequency .* h_Wnorm_frequency;
filtered_noisy_x_Wnorm = VW * filtered_noisy_x_Wnorm_frequency;
filtered_noisy_x_L_frequency = noisy_x_L_frequency .* h_L_frequency;
filtered_noisy_x_L = G.U * filtered_noisy_x_L_frequency;
figure();
gsp_plot_signal(G, x);
title('x');
figure();
gsp_plot_signal(G, noisy_x);
title('noisy x');
figure();
gsp_plot_signal(G, filtered_noisy_x_Wnorm);
title('filtered noisy x (Wnorm)');
figure();
gsp_plot_signal(G, filtered_noisy_x_L);
title('filtered noisy x (L)');
%% 7)
noise_P_Wnorm = 1 / length(x) * (filtered_noisy_x_Wnorm - x)' * (filtered_noisy_x_Wnorm - x);
SNR_Wnorm = 10 * log10(signal_P / noise_P_Wnorm);
noise_P_L = 1 / length(x) * (filtered_noisy_x_L - x)' * (filtered_noisy_x_L - x);
SNR_L = 10 * log10(signal_P / noise_P_L);
%% 8)
lambda_Wnorm = zeros(length(diag(DWnorm)));
for i = 1:length(diag(DWnorm))
    lambda_Wnorm(:, i) = diag(DWnorm) .^ (i - 1);
end
h_FIR_Wnorm = lambda_Wnorm \ h_Wnorm_frequency;
lambda_L = zeros(length(G.e));
for i = 1:length(G.e)
    lambda_L(:, i) = G.e .^ (i - 1);
end
h_FIR_L = lambda_L\ h_L_frequency;
%% 9)
h_FIR_Wnorm_approximate = (lambda_Wnorm(:, 1:3)' * lambda_Wnorm(:, 1:3)) \ (lambda_Wnorm(:, 1:3)' * h_Wnorm_frequency);
h_FIR_L_approximate = (lambda_L(:, 1:3)' * lambda_L(:, 1:3)) \ (lambda_L(:, 1:3)' * h_L_frequency);
%% 10)
figure();
stem(diag(DWnorm), h_Wnorm_frequency)
xlim([-1 1])
ylim([0 1])
title('frequency response of h (Wnorm)');
h_Wnorm_frequency_approximate = lambda_Wnorm(:, 1:3) * h_FIR_Wnorm_approximate;
figure();
stem(diag(DWnorm), h_Wnorm_frequency_approximate)
xlim([-1 1])
ylim([0 1])
title('frequency response of approximate h (Wnorm)');
figure();
stem(G.e, h_L_frequency)
xlim([0 8])
ylim([0 1])
title('frequency response of h (L)');
h_L_frequency_approximate = lambda_L(:, 1:3) * h_FIR_L_approximate;
figure();
stem(G.e, h_L_frequency_approximate)
xlim([0 8])
ylim([0 1])
title('frequency response of approximate h (L)');
%% 11)
filtered_noisy_x_Wnorm_frequency_approximate = noisy_x_Wnorm_frequency .* h_Wnorm_frequency_approximate;
filtered_noisy_x_Wnorm_approximate = VW * filtered_noisy_x_Wnorm_frequency_approximate;
filtered_noisy_x_L_frequency_approximate = noisy_x_L_frequency .* h_L_frequency_approximate;
filtered_noisy_x_L_approximate = G.U * filtered_noisy_x_L_frequency_approximate;
figure();
gsp_plot_signal(G, x);
title('x');
figure();
gsp_plot_signal(G, noisy_x);
title('noisy x');
figure();
gsp_plot_signal(G, filtered_noisy_x_Wnorm_approximate);
title('approximately filtered noisy x (Wnorm)');
figure();
gsp_plot_signal(G, filtered_noisy_x_L_approximate);
title('approximately filtered noisy x (L)');
noise_P_Wnorm_approximate = 1 / length(x) * (filtered_noisy_x_Wnorm_approximate - x)' * (filtered_noisy_x_Wnorm_approximate - x);
SNR_Wnorm_approximate = 10 * log10(signal_P / noise_P_Wnorm_approximate);
noise_P_L_approximate = 1 / length(x) * (filtered_noisy_x_L_approximate - x)' * (filtered_noisy_x_L_approximate - x);
SNR_L_approximate = 10 * log10(signal_P / noise_P_L_approximate);
%% 12)
h_FIR_Wnorm_approximate_2 = (lambda_Wnorm(:, 1:2)' * lambda_Wnorm(:, 1:2)) \ (lambda_Wnorm(:, 1:2)' * h_Wnorm_frequency);
h_FIR_L_approximate_2 = (lambda_L(:, 1:2)' * lambda_L(:, 1:2)) \ (lambda_L(:, 1:2)' * h_L_frequency);
h_Wnorm_frequency_approximate_2 = lambda_Wnorm(:, 1:2) * h_FIR_Wnorm_approximate_2;
figure();
stem(diag(DWnorm), h_Wnorm_frequency_approximate_2)
xlim([-1 1])
ylim([0 1])
title('frequency response of approximate h with 2 coefficients (Wnorm)');
h_L_frequency_approximate_2 = lambda_L(:, 1:2) * h_FIR_L_approximate_2;
figure();
stem(G.e, h_L_frequency_approximate_2)
xlim([0 8])
ylim([0 1])
title('frequency response of approximate h with 2 coefficients (L)');
filtered_noisy_x_Wnorm_frequency_approximate_2 = noisy_x_Wnorm_frequency .* h_Wnorm_frequency_approximate_2;
filtered_noisy_x_Wnorm_approximate_2 = VW * filtered_noisy_x_Wnorm_frequency_approximate_2;
filtered_noisy_x_L_frequency_approximate_2 = noisy_x_L_frequency .* h_L_frequency_approximate_2;
filtered_noisy_x_L_approximate_2 = G.U * filtered_noisy_x_L_frequency_approximate_2;
figure();
gsp_plot_signal(G, filtered_noisy_x_Wnorm_approximate_2);
title('approximately filtered noisy x with 2 coefficients (Wnorm)');
figure();
gsp_plot_signal(G, filtered_noisy_x_L_approximate_2);
title('approximately filtered noisy x with 2 coefficients (L)');
noise_P_Wnorm_approximate_2 = 1 / length(x) * (filtered_noisy_x_Wnorm_approximate_2 - x)' * (filtered_noisy_x_Wnorm_approximate_2 - x);
SNR_Wnorm_approximate_2 = 10 * log10(signal_P / noise_P_Wnorm_approximate_2);
noise_P_L_approximate_2 = 1 / length(x) * (filtered_noisy_x_L_approximate_2 - x)' * (filtered_noisy_x_L_approximate_2 - x);
SNR_L_approximate_2 = 10 * log10(signal_P / noise_P_L_approximate_2);
h_FIR_Wnorm_approximate_5 = (lambda_Wnorm(:, 1:5)' * lambda_Wnorm(:, 1:5)) \ (lambda_Wnorm(:, 1:5)' * h_Wnorm_frequency);
h_FIR_L_approximate_5 = (lambda_L(:, 1:5)' * lambda_L(:, 1:5)) \ (lambda_L(:, 1:5)' * h_L_frequency);
h_Wnorm_frequency_approximate_5 = lambda_Wnorm(:, 1:5) * h_FIR_Wnorm_approximate_5;
figure();
stem(diag(DWnorm), h_Wnorm_frequency_approximate_5)
xlim([-1 1])
ylim([0 1])
title('frequency response of approximate h with 5 coefficients (Wnorm)');
h_L_frequency_approximate_5 = lambda_L(:, 1:5) * h_FIR_L_approximate_5;
figure();
stem(G.e, h_L_frequency_approximate_5)
xlim([0 8])
ylim([0 1])
title('frequency response of approximate h with 5 coefficients (L)');
filtered_noisy_x_Wnorm_frequency_approximate_5 = noisy_x_Wnorm_frequency .* h_Wnorm_frequency_approximate_5;
filtered_noisy_x_Wnorm_approximate_5 = VW * filtered_noisy_x_Wnorm_frequency_approximate_5;
filtered_noisy_x_L_frequency_approximate_5 = noisy_x_L_frequency .* h_L_frequency_approximate_5;
filtered_noisy_x_L_approximate_5 = G.U * filtered_noisy_x_L_frequency_approximate_5;
figure();
gsp_plot_signal(G, filtered_noisy_x_Wnorm_approximate_5);
title('approximately filtered noisy x with 5 coefficients (Wnorm)');
figure();
gsp_plot_signal(G, filtered_noisy_x_L_approximate_5);
title('approximately filtered noisy x with 5 coefficients (L)');
noise_P_Wnorm_approximate_5 = 1 / length(x) * (filtered_noisy_x_Wnorm_approximate_5 - x)' * (filtered_noisy_x_Wnorm_approximate_5 - x);
SNR_Wnorm_approximate_5 = 10 * log10(signal_P / noise_P_Wnorm_approximate_5);
noise_P_L_approximate_5 = 1 / length(x) * (filtered_noisy_x_L_approximate_5 - x)' * (filtered_noisy_x_L_approximate_5 - x);
SNR_L_approximate_5 = 10 * log10(signal_P / noise_P_L_approximate_5);
%% b)
%% 1)
A = [0 1 1 1 1 1 1 1; 1 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0;
     1 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0];
W = A;
C = [0 0; cos(0 * pi / 7) sin(0 * pi / 7); cos(2 * pi / 7) sin(2 * pi / 7); cos(4 * pi / 7) sin(4 * pi / 7);
    cos(6 * pi / 7) sin(6 * pi / 7); cos(8 * pi / 7) sin(8 * pi / 7); cos(10 * pi / 7) sin(10 * pi / 7); cos(12 * pi / 7) sin(12 * pi / 7)];
star_G = gsp_graph(W, C);
star_G.plotting.edge_color = 'b';
star_G.plotting.vertex_color = 'k';
signal = zeros(8, 1);
signal(1) = 100;
figure();
gsp_plot_signal(star_G, signal);
title('Star Graph');
S = eye(8);
[VW,DW] = eig(W);
Wnorm = 1 / max(abs(DW), [], 'all') * W;
DWnorm = DW / max(abs(DW), [], 'all');
star_x_Wnorm_frequency = VW \ signal;
star_S_Wnorm_frequency = VW \ S;
star_shifted_signals_Wnorm_frequency = star_S_Wnorm_frequency .* star_x_Wnorm_frequency;
star_shifted_signals_Wnorm = VW * star_shifted_signals_Wnorm_frequency;
star_G = gsp_compute_fourier_basis(star_G);
star_x_L_frequency = star_G.U \ signal;
star_S_L_frequency = star_G.U \ S;
star_shifted_signals_L_frequency = star_S_L_frequency .* star_x_L_frequency;
star_shifted_signals_L = star_G.U * star_shifted_signals_L_frequency;
%% 2)
A = [0 1 0 0 0 0 0 1; 1 0 1 0 0 0 0 0; 0 1 0 1 0 0 0 0; 0 0 1 0 1 0 0 0;
    0 0 0 1 0 1 0 0; 0 0 0 0 1 0 1 0; 0 0 0 0 0 1 0 1; 1 0 0 0 0 0 1 0];
W = A;
C = [cos(0 * pi / 8) sin(0 * pi / 8); cos(2 * pi / 8) sin(2 * pi / 8); cos(4 * pi / 8) sin(4 * pi / 8); cos(6 * pi / 8) sin(6 * pi / 8);
    cos(8 * pi / 8) sin(8 * pi / 8); cos(10 * pi / 8) sin(10 * pi / 8); cos(12 * pi / 8) sin(12 * pi / 8); cos(14 * pi / 8) sin(14 * pi / 8)];
cycle_G = gsp_graph(W, C);
cycle_G.plotting.edge_color = 'b';
cycle_G.plotting.vertex_color = 'k';
signal = zeros(8, 1);
signal(1) = 100;
figure();
gsp_plot_signal(cycle_G, signal);
title('Cycle Graph');
S = eye(8);
[VW,DW] = eig(W);
Wnorm = 1 / max(abs(DW), [], 'all') * W;
DWnorm = DW / max(abs(DW), [], 'all');
cycle_x_Wnorm_frequency = VW \ signal;
cycle_S_Wnorm_frequency = VW \ S;
cycle_shifted_signals_Wnorm_frequency = cycle_S_Wnorm_frequency .* cycle_x_Wnorm_frequency;
cycle_shifted_signals_Wnorm = VW * cycle_shifted_signals_Wnorm_frequency;
cycle_G = gsp_compute_fourier_basis(cycle_G);
cycle_x_L_frequency = cycle_G.U \ signal;
cycle_S_L_frequency = cycle_G.U \ S;
cycle_shifted_signals_L_frequency = cycle_S_L_frequency .* cycle_x_L_frequency;
cycle_shifted_signals_L = cycle_G.U * cycle_shifted_signals_L_frequency;
%% 3)
A = [0 1 1 1 1 1 1 1; 1 0 1 1 1 1 1 1; 1 1 0 1 1 1 1 1; 1 1 1 0 1 1 1 1;
    1 1 1 1 0 1 1 1; 1 1 1 1 1 0 1 1; 1 1 1 1 1 1 0 1; 1 1 1 1 1 1 1 0];
W = A;
C = [cos(0 * pi / 8) sin(0 * pi / 8); cos(2 * pi / 8) sin(2 * pi / 8); cos(4 * pi / 8) sin(4 * pi / 8); cos(6 * pi / 8) sin(6 * pi / 8);
    cos(8 * pi / 8) sin(8 * pi / 8); cos(10 * pi / 8) sin(10 * pi / 8); cos(12 * pi / 8) sin(12 * pi / 8); cos(14 * pi / 8) sin(14 * pi / 8)];
complete_G = gsp_graph(W, C);
complete_G.plotting.edge_color = 'b';
complete_G.plotting.vertex_color = 'k';
signal = zeros(8, 1);
signal(1) = 100;
figure();
gsp_plot_signal(complete_G, signal);
title('Complete Graph');
S = eye(8);
[VW,DW] = eig(W);
Wnorm = 1 / max(abs(DW), [], 'all') * W;
DWnorm = DW / max(abs(DW), [], 'all');
complete_x_Wnorm_frequency = VW \ signal;
complete_S_Wnorm_frequency = VW \ S;
complete_shifted_signals_Wnorm_frequency = complete_S_Wnorm_frequency .* complete_x_Wnorm_frequency;
complete_shifted_signals_Wnorm = VW * complete_shifted_signals_Wnorm_frequency;
complete_G = gsp_compute_fourier_basis(complete_G);
complete_x_L_frequency = complete_G.U \ signal;
complete_S_L_frequency = complete_G.U \ S;
complete_shifted_signals_L_frequency = complete_S_L_frequency .* complete_x_L_frequency;
complete_shifted_signals_L = complete_G.U * complete_shifted_signals_L_frequency;