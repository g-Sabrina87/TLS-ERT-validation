%% 10_3_make_figure3_block_size_sensitivity.m
% Figure 3 of the paper: differences in hold-out AUC and IoU from the unconstrained
% model U as a function of the horizontal size of the spatial blocks (error bars:
% 95% bootstrap intervals), for H (panels a, b) and C (panels c, d).
% Values: Table 3 of the paper (medians and 95% limits of the paired spatial block
% bootstrap; blocks of 3, 5, 7 and 10 m; the 15 and 20 m blocks added at Chiusazza
% are in Table 3 only). Columns: [median, low, high].
% Output: Figure_block_size_sensitivity.png

clear; clc; close all;

blocks = [3 5 7 10];
%            dAUC H-U                 dIoU H-U                  dAUC C-U                dIoU C-U
M.dAUC_H = [ 0.000 -0.006  0.006;  0.000 -0.008  0.007;  0.000 -0.006  0.006; -0.001 -0.005  0.008];
M.dIoU_H = [-0.066 -0.125 -0.016; -0.064 -0.134 -0.010; -0.066 -0.126 -0.013; -0.064 -0.160  0.000];
M.dAUC_C = [ 0.192  0.123  0.269;  0.193  0.096  0.308;  0.189  0.096  0.315;  0.187  0.060  0.361];
M.dIoU_C = [ 0.155  0.070  0.257;  0.155  0.047  0.287;  0.152  0.059  0.251;  0.155  0.036  0.347];
C.dAUC_H = [ 0.007  0.003  0.010;  0.007  0.002  0.011;  0.007  0.001  0.012;  0.006  0.000  0.013];
C.dIoU_H = [ 0.009 -0.003  0.024;  0.009 -0.006  0.026;  0.010 -0.008  0.029;  0.009 -0.011  0.035];
C.dAUC_C = [ 0.210  0.178  0.245;  0.210  0.165  0.263;  0.208  0.152  0.271;  0.209  0.148  0.279];
C.dIoU_C = [ 0.542  0.500  0.581;  0.541  0.485  0.595;  0.542  0.480  0.596;  0.541  0.465  0.608];

cMicio = [0.75 0.25 0.00]; cChiu = [0.17 0.45 0.62];
panels = {'dAUC_H', '(a) H vs U: \DeltaAUC'; 'dIoU_H', '(b) H vs U: \DeltaIoU'; ...
          'dAUC_C', '(c) C vs U: \DeltaAUC'; 'dIoU_C', '(d) C vs U: \DeltaIoU'};
dx = 0.15 * blocks;          % small horizontal offset so that the sites do not overlap

figure('Color', 'w', 'Position', [60 40 1250 900]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
for k = 1:4
    nexttile; hold on; box on;
    A = C.(panels{k,1}); B = M.(panels{k,1});
    errorbar(blocks - dx/2, A(:,1), A(:,1)-A(:,2), A(:,3)-A(:,1), '-o', 'Color', cChiu, ...
             'MarkerFaceColor', cChiu, 'LineWidth', 1.4, 'CapSize', 6, 'DisplayName', 'Chiusazza');
    errorbar(blocks + dx/2, B(:,1), B(:,1)-B(:,2), B(:,3)-B(:,1), '-s', 'Color', cMicio, ...
             'MarkerFaceColor', cMicio, 'LineWidth', 1.4, 'CapSize', 6, 'DisplayName', 'Micio Conti');
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    set(gca, 'XTick', blocks); xlim([2 11.5]);
    xlabel('Block size, horizontal (m)');
    if k == 1 || k == 3, ylabel('Difference (model - U)'); end
    title(panels{k,2}, 'FontWeight', 'normal');
    legend('Location', 'best');
end
exportgraphics(gcf, 'Figure_block_size_sensitivity.png', 'Resolution', 300);
