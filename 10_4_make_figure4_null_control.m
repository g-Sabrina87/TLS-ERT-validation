%% 10_4_make_figure4_null_control.m
% Figure 4 of the paper: AUC difference from U on the main evaluation set for the
% null constraints (dots), their mean (solid line) and the correct constraint H
% (dashed line), at Micio Conti (a) and Chiusazza (b).
% Inputs: the CSV files written by null_diagnostics_buffer.m (columns: model, AUC_rank,
% IoU, nAboveThr, AUC_grid; rows U, H, NULL01, ...), for the main evaluation set:
%   null_buffer_Micio_Conti_excl1.5_near3.csv   (exclusion 1.5 m, distance 3 m)
%   null_buffer_Chiusazza_excl1_near3.csv       (exclusion 1.0 m, distance 3 m)
% Output: Figure_null_control.png

clear; clc; close all;

files  = {'null_buffer_Micio_Conti_excl1.5_near3.csv', 'null_buffer_Chiusazza_excl1_near3.csv'};
titles = {'(a) Micio Conti', '(b) Chiusazza'};
cols   = {[0.75 0.25 0.00], [0.17 0.45 0.62]};
rng(1);                                     % horizontal jitter of the dots (display only)

figure('Color', 'w', 'Position', [60 60 1250 560]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    T = readtable(files{s});
    names = string(T.model);
    aucU = T.AUC_rank(names == "U");
    aucH = T.AUC_rank(names == "H");
    isNull = startsWith(names, "NULL");
    d = T.AUC_rank(isNull) - aucU;          % AUC difference of each null from U
    dH = aucH - aucU;
    nNull = numel(d);
    nGE = sum(T.AUC_rank(isNull) >= aucH);  % nulls with AUC >= H
    fprintf('%s: %d nulls, mean difference %+.5f, H %+.5f, nulls >= H: %d\n', titles{s}, nNull, mean(d), dH, nGE);

    nexttile; hold on; box on;
    xj = 1 + 0.35*(rand(nNull,1) - 0.5);
    plot(xj, d, 'o', 'MarkerFaceColor', cols{s}, 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
         'DisplayName', sprintf('%d null constraints', nNull));
    plot([0.55 1.45], [mean(d) mean(d)], '-', 'Color', cols{s}, 'LineWidth', 2, 'DisplayName', 'Null mean');
    plot([0.55 1.45], [dH dH], 'k--', 'LineWidth', 2.2, 'DisplayName', 'Correct constraint (H)');
    yline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); text(1.47, 0, 'U', 'Color', [0.5 0.5 0.5], 'VerticalAlignment', 'middle');
    xlim([0.5 1.5]); set(gca, 'XTick', []);
    ylabel('AUC difference from U (hold-out, main evaluation set)');
    title(titles{s}, 'FontWeight', 'normal');
    text(1.45, min(ylim) + 0.04*diff(ylim), sprintf('%d of %d nulls \\geq H', nGE, nNull), ...
         'HorizontalAlignment', 'right');
    legend('Location', 'southoutside', 'Orientation', 'horizontal');
end
exportgraphics(gcf, 'Figure_null_control.png', 'Resolution', 300);
