%% 1)
gsp_start;
%% 2)
W1 = [0, 0.7, 1.1, 2.3; 0.7, 0, 0, 0; 1.1, 0, 0, 0; 2.3, 0, 0, 0];
C1 = [0, 0; -1, -sqrt(3); 1, -sqrt(3); 0, 2];
G1 = gsp_graph(W1, C1);
G1.plotting.edge_color = 'b';
G1.plotting.vertex_color = 'k';
figure();
gsp_plot_graph(G1);
W2 = [0, 1.6, 2.4; 1.6, 0, 0.8; 2.4, 0.8, 0];
C2 = [0, 1; -1, -1; 1, -1];
G2 = gsp_graph(W2, C2);
G2.plotting.edge_color = 'b';
G2.plotting.vertex_color = 'k';
figure();
gsp_plot_graph(G2);
%% 3)
C = zeros(size(C1,1)*size(C2,1), 2);
for i=1:size(C1,1)
    C((i-1)*size(C2,1)+1:i*size(C2,1), 1) = C2(:, 1) + C1(i, 1)/5;
    C((i-1)*size(C2,1)+1:i*size(C2,1), 2) = C2(:, 2) + C1(i, 2)/5;
end
param.verbose = 1;
param.rule = 'kronecker';
Gt = gsp_graph_product(G1, G2, param);
Gt.coords = C;
Gt_Weights = Gt.W;
Gt_Adjacency = Gt.A;
Gt.plotting.edge_color = 'b';
Gt.plotting.vertex_color = 'k';
figure();
gsp_plot_graph(Gt);
param.verbose = 1;
param.rule = 'cartesian';
Gs = gsp_graph_product(G1, G2, param);
Gs.coords = C;
Gs_Weights = Gs.W;
Gs_Adjacency = Gs.A;
Gs.plotting.edge_color = 'b';
Gs.plotting.vertex_color = 'k';
figure();
gsp_plot_graph(Gs);
%% 4)
param.verbose = 1;
param.rule = 'cartesian';
myG = gsp_graph_product(G1, G2, param);
myG.coords = C;
signal = (rand(myG.N, 1) - 0.5) * 20;
figure();
gsp_plot_signal(myG, signal);
%% 5)
myG = gsp_compute_fourier_basis(myG);
figure();
stem(1:myG.N, myG.e);
%% 6)
figure();
gsp_plot_signal(myG, myG.U(:, 1));
figure();
gsp_plot_signal(myG, myG.U(:, 2));
figure();
gsp_plot_signal(myG, myG.U(:, 11));
figure();
gsp_plot_signal(myG, myG.U(:, 12));
%% 7)
GL = gsp_logo;
signal = zeros(GL.N, 1);
for i=1:GL.N
    if GL.coords(i, 1) <= 200
        signal(i, 1) = -1;
    elseif GL.coords(i, 1) <= 400
        signal(i, 1) = 1;
    else
        signal(i, 1) = -0.5;
    end
end
figure();
gsp_plot_signal(GL, signal);
%% 8)
GL = gsp_compute_fourier_basis(GL);
%% 9)
signal_2D = GL.U(:, 2:3);
%% 10)
figure();
scatter(signal_2D(:, 1), signal_2D(:, 2), 1);
%% 11)
idx2 = kmeans(signal_2D, 3);
%% 12)
figure();
gsp_plot_signal(GL, idx2);
%% 13)
signal_3D = GL.U(:, 2:4);
figure();
scatter3(signal_3D(:, 1), signal_3D(:, 2), signal_3D(:, 3), 1);
idx3 = kmeans(signal_3D, 3);
figure();
gsp_plot_signal(GL, idx3);