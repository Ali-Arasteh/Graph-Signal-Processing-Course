gsp_start;
%% b)
%% 10)
[~,~,data] = xlsread('Data_city.csv');
temp = [];
for i=1:length(data)
    if strcmp(data{i,7},'primary') ||  strcmp(data{i,7},'admin')
        temp = [temp, i];
    end
end
data = data(temp, [2, 3, 10]);
distances = zeros(length(data), length(data));
for i=1:length(data)
    for j=i+1:length(data)
       distances(i, j) = exp(-(getDistance(data{i, 1}, data{i, 2}, data{j, 1}, data{j, 2})/500)^2);
       if distances(i, j) <= 0.2
           distances(i, j) = 0;
       end
       distances(j, i) = distances(i, j);
    end
end
radius = 6371;
locations = zeros(length(data), 2);
locations(:, 1) = radius*cell2mat(data(:, 2));
locations(:, 2) = radius*cell2mat(data(:, 1));
G = gsp_graph(distances, locations);
G.plotting.edge_color = 'b';
G.plotting.vertex_color = 'k';
figure();
gsp_plot_graph(G);
title('G');
%% 11)
tempreture = cell2mat(data(:, 3));
figure();
gsp_plot_signal(G, tempreture);
title('tempreture');
[indices_cell_1, sk_reduced_G_cell_1, downsampled_filtered_x_cell_1, interpolation_error_cell_1] = MyPyramidAnalysis(G, tempreture, 3, 300, true, 0);
figure();
gsp_plot_signal(sk_reduced_G_cell_1{2}, downsampled_filtered_x_cell_1{2});
title('tempreture in 1 level reduction - way 1');
figure();
gsp_plot_signal(sk_reduced_G_cell_1{1}, interpolation_error_cell_1{1});
title('interpolation error in 1 level reducrion - way 1');
figure();
gsp_plot_signal(sk_reduced_G_cell_1{3}, downsampled_filtered_x_cell_1{3});
title('tempreture in 2 level reduction - way 1');
figure();
gsp_plot_signal(sk_reduced_G_cell_1{2}, interpolation_error_cell_1{2});
title('interpolation error in 2 level reducrion - way 1');
figure();
gsp_plot_signal(sk_reduced_G_cell_1{4}, downsampled_filtered_x_cell_1{4});
title('tempreture in 3 level reduction - way 1');
figure();
gsp_plot_signal(sk_reduced_G_cell_1{3}, interpolation_error_cell_1{3});
title('interpolation error in 3 level reducrion - way 1');
total_interpolation_error_1 = zeros(3, 1);
for i=1:3
    total_interpolation_error_1(i) = norm(interpolation_error_cell_1{i});
end
average_interpolation_error_1 = sum(total_interpolation_error_1, 'all')/(length(interpolation_error_cell_1{1})+length(interpolation_error_cell_1{2})+length(interpolation_error_cell_1{3}));
[indices_cell_2, sk_reduced_G_cell_2, downsampled_filtered_x_cell_2, interpolation_error_cell_2] = MyPyramidAnalysis(G, tempreture, 3, 300, false, 0.1);
figure();
gsp_plot_signal(sk_reduced_G_cell_2{2}, downsampled_filtered_x_cell_2{2});
title('tempreture in 1 level reduction - way 2');
figure();
gsp_plot_signal(sk_reduced_G_cell_2{1}, interpolation_error_cell_2{1});
title('interpolation error in 1 level reducrion - way 2');
figure();
gsp_plot_signal(sk_reduced_G_cell_2{3}, downsampled_filtered_x_cell_2{3});
title('tempreture in 2 level reduction - way 2');
figure();
gsp_plot_signal(sk_reduced_G_cell_2{2}, interpolation_error_cell_2{2});
title('interpolation error in 2 level reducrion - way 2');
figure();
gsp_plot_signal(sk_reduced_G_cell_2{4}, downsampled_filtered_x_cell_2{4});
title('tempreture in 3 level reduction - way 2');
figure();
gsp_plot_signal(sk_reduced_G_cell_2{3}, interpolation_error_cell_2{3});
title('interpolation error in 3 level reducrion - way 2');
total_interpolation_error_2 = zeros(3, 1);
for i=1:3
    total_interpolation_error_2(i) = norm(interpolation_error_cell_2{i});
end
average_interpolation_error_2 = sum(total_interpolation_error_2, 'all')/(length(interpolation_error_cell_2{1})+length(interpolation_error_cell_2{2})+length(interpolation_error_cell_2{3}));
%% 12)
figure();
gsp_plot_signal(G, tempreture);
title('tempreture');
tempreture_1 = MyPyramidSynthesis(downsampled_filtered_x_cell_1{4}, indices_cell_1, sk_reduced_G_cell_1, interpolation_error_cell_1, 3, true, 0);
figure();
gsp_plot_signal(G, tempreture_1);
title('recovered tempreture - way 1');
tempreture_2 = MyPyramidSynthesis(downsampled_filtered_x_cell_2{4}, indices_cell_2, sk_reduced_G_cell_2, interpolation_error_cell_2, 3, false, 0.1);
figure();
gsp_plot_signal(G, tempreture_2);
title('recovered tempreture - way 2');
%% 13
figure();
gsp_plot_signal(G, tempreture);
title('tempreture');
SNR = 5;
signal_P_1 = 1 / length(downsampled_filtered_x_cell_1{4}) * (downsampled_filtered_x_cell_1{4}') * downsampled_filtered_x_cell_1{4};
noise_variance_1 = signal_P_1 / 10^(SNR / 10);
noisy_downsampled_filtered_x_cell_1 = downsampled_filtered_x_cell_1{4} + randn(length(downsampled_filtered_x_cell_1{4}), 1) * sqrt(noise_variance_1);
noisy_tempreture_1 = MyPyramidSynthesis(noisy_downsampled_filtered_x_cell_1, indices_cell_1, sk_reduced_G_cell_1, interpolation_error_cell_1, 3, true, 0);
figure();
gsp_plot_signal(G, noisy_tempreture_1);
title('recovered tempreture - way 1');
MSE_1 = mse(tempreture-noisy_tempreture_1);
signal_P_2 = 1 / length(downsampled_filtered_x_cell_2{4}) * (downsampled_filtered_x_cell_2{4}') * downsampled_filtered_x_cell_2{4};
noise_variance_2 = signal_P_2 / 10^(SNR / 10);
noisy_downsampled_filtered_x_cell_2 = downsampled_filtered_x_cell_2{4} + randn(length(downsampled_filtered_x_cell_2{4}), 1) * sqrt(noise_variance_2);
noisy_tempreture_2 = MyPyramidSynthesis(noisy_downsampled_filtered_x_cell_2, indices_cell_2, sk_reduced_G_cell_2, interpolation_error_cell_2, 3, false, 0.1);
figure();
gsp_plot_signal(G, noisy_tempreture_2);
title('recovered tempreture - way 2');
MSE_2 = mse(tempreture-noisy_tempreture_2);
%% c)
%% 14)
G = gsp_compute_fourier_basis(G);
c = 0.01;
temp = G.U*diag(exp(-c*G.e))*pinv(G.U);
tempreture_0 = tempreture;
tempreture_1 = temp*tempreture_0;
figure();
gsp_plot_signal(G, tempreture_1);
title('tempreture on graph - day 1');
figure();
stem(1:length(tempreture_0), tempreture_1);
title('tempreture - day 1');
tempreture_2 = temp*tempreture_1;
figure();
gsp_plot_signal(G, tempreture_2);
title('tempreture on graph - day 2');
figure();
stem(1:length(tempreture_0), tempreture_2);
title('tempreture - day 2');
tempreture_3 = temp*tempreture_2;
figure();
gsp_plot_signal(G, tempreture_3);
title('tempreture on graph - day 3');
figure();
stem(1:length(tempreture_0), tempreture_3);
title('tempreture - day 3');
tempreture_4 = temp*tempreture_3;
figure();
gsp_plot_signal(G, tempreture_4);
title('tempreture on graph - day 4');
figure();
stem(1:length(tempreture_0), tempreture_4);
title('tempreture - day 4');
tempreture_5 = temp*tempreture_4;
figure();
gsp_plot_signal(G, tempreture_5);
title('tempreture on graph - day 5');
figure();
stem(1:length(tempreture_0), tempreture_5);
title('tempreture - day 5');
lambda_inf = zeros(length(tempreture_0), 1);
lambda_inf(1) = 1;
tempreture_inf = G.U*diag(lambda_inf)*pinv(G.U)*tempreture_0;
figure();
gsp_plot_signal(G, tempreture_inf);
title('tempreture on graph - day infinity');
figure();
stem(1:length(tempreture_0), tempreture_inf);
title('tempreture - day infinity');
%% a) functions
%% 1) 
function indices = MyVertexSelection(G)
G = gsp_compute_fourier_basis(G);
umax = G.U(:, end);
indices = [];
for i=1:length(umax)
    if (umax(i) >= 0)
        indices = [indices, i];
    end
end
end
%% 2)
function sk_reduced_G = MySKReduction(G, indices, Q)
indices_c = setdiff(1:G.N, indices);
kron_reduced_L = G.L(indices, indices)-G.L(indices, indices_c)*pinv(G.L(indices_c, indices_c))*G.L(indices_c, indices);
kron_reduced_W = diag(diag(kron_reduced_L))-kron_reduced_L;
dRG = zeros(length(indices), length(indices));
for i=1:length(indices)
    for j=1:length(indices)
        delta = zeros(length(indices), 1);
        delta(i) = delta(i)+1;
        delta(j) = delta(j)-1;
        dRG(i, j) = delta'*pinv(kron_reduced_L)*delta;
    end
end
P = dRG.*kron_reduced_W/sum(dRG.*kron_reduced_W, 'all');
sk_reduced_W = zeros(length(indices), length(indices));
for i=1:Q
    temp = randsrc(1,1,[1:length(indices)^2;reshape(P, 1, [])]);
    m = floor((temp-1)/length(indices))+1;
    n = rem((temp-1), length(indices))+1;
    sk_reduced_W(m, n) = sk_reduced_W(m, n)+kron_reduced_W(m, n)/(Q*2*P(m, n));
    sk_reduced_W(n, m) = sk_reduced_W(m, n);
end
sk_reduced_G = gsp_graph(sk_reduced_W, G.coords(indices, :));
end
%% 3)
function filtered_x = MyHfilter(G, x)
G = gsp_compute_fourier_basis(G);
h_frequency = 1./(1+2*G.e);
x_frequency = G.U\x;
filtered_x_frequency = h_frequency.*x_frequency;
filtered_x = G.U*filtered_x_frequency;
end
%% 4)
function downsampled_x = MyDS(x, indices)
downsampled_x = x(indices);
end
%% 5)
function interpolated_x = MyInterpolate(G, indices, downsampled_x, flag, epsilon)
G = gsp_compute_fourier_basis(G);
if flag
    interpolated_x = G.U*pinv(G.U(indices, :))*downsampled_x;
%     temp1 = pinv(G.U(indices, :))*downsampled_x;
%     temp2 = temp1.^2/sum(temp1.^2, 'all');
%     i = 0;
%     s = 0;
%     m = zeros(length(temp1), 1);
%     while s < 0.75
%         i = i+1;
%         s = s+temp2(i);
%         m(i) = 1;
%     end
%     temp1 = temp1.*m;
%     interpolated_x = G.U*temp1;
else
    regularized_L = G.L+epsilon*eye(G.N);
    Phi = pinv(regularized_L);
    interpolated_x = Phi(:, indices)*pinv(Phi(indices, indices))*downsampled_x;
end
end
%% 6)
function [indices, sk_reduced_G, downsampled_filtered_x, interpolation_error] = MyAnalysis(G, x, Q, flag, epsilon)
indices = MyVertexSelection(G);
sk_reduced_G = MySKReduction(G, indices, Q);
filtered_x = MyHfilter(G, x);
downsampled_filtered_x = MyDS(filtered_x, indices);
interpolated_filtered_x = MyInterpolate(G, indices, downsampled_filtered_x, flag, epsilon);
interpolation_error = x-interpolated_filtered_x;
end
%% 7)
function [indices_cell, sk_reduced_G_cell, downsampled_filtered_x_cell, interpolation_error_cell] = MyPyramidAnalysis(G, x, N, Q, flag, epsilon)
indices_cell = cell(N,1);
sk_reduced_G_cell = cell(N+1, 1);
downsampled_filtered_x_cell = cell(N+1, 1);
interpolation_error_cell = cell(N, 1);
sk_reduced_G_cell{1} = G;
downsampled_filtered_x_cell{1} = x;
for i=1:N
    [indices_cell{i}, sk_reduced_G_cell{i+1}, downsampled_filtered_x_cell{i+1}, interpolation_error_cell{i}] = MyAnalysis(G, x, floor(Q/2^(i-1)), flag, epsilon);
    G = sk_reduced_G_cell{i+1};
    x = downsampled_filtered_x_cell{i+1};
end
end
%% 8)
function x = MySynthesis(downsampled_filtered_x, indices, sk_reduced_G, interpolation_error, flag, epsilon)
x = MyInterpolate(sk_reduced_G, indices, downsampled_filtered_x, flag, epsilon)+interpolation_error;
end
%% 9)
function x = MyPyramidSynthesis(downsampled_filtered_x, indices_cell, sk_reduced_G_cell, interpolation_error_cell, N, flag, epsilon)
for i=1:N
    x = MySynthesis(downsampled_filtered_x, indices_cell{N-i+1}, sk_reduced_G_cell{N-i+1}, interpolation_error_cell{N-i+1}, flag, epsilon);
    downsampled_filtered_x = x;
end
end