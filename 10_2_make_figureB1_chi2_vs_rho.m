%% 10_2_make_figure2_chi2_vs_rho.m
% Figure 2 of the paper: reduced chi-squared as a function of the imposed
% resistivity, for Micio Conti (a) and Chiusazza (b).
% The values are those of Table 1 of the paper (3 decimals). If you have the
% full-precision values from the ERTLab+ logs (final Data Resid / number of
% quadrupoles: 1938 at Micio Conti, 3535 at Chiusazza), replace them below.
% At Chiusazza a saturation curve chi2 = a + b*exp(-rho/tau) is fitted to the five
% constrained values (fminsearch, base MATLAB); at Micio Conti the four values
% span only 0.004 and no curve is fitted.
% Output: Figure_chi2_saturation.png

clear; clc; close all;

cMicio = [0.75 0.25 0.00]; cChiu = [0.17 0.45 0.62];

% ---- data (Table 1) ----
rhoM = [17000 50000 75000 100000];       chiM = [0.601 0.599 0.598 0.597];   freeM = 0.615;
rhoC = [4000 10000 20000 50000 100000];  chiC = [0.367 0.347 0.339 0.333 0.324]; freeC = 0.419;

% ---- saturation fit at Chiusazza ----
model = @(p, r) p(1) + p(2).*exp(-r./p(3));
cost  = @(p) sum((model(p, rhoC) - chiC).^2);
p0    = [min(chiC), chiC(1)-min(chiC), mean(rhoC)/3];
opts  = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'TolFun', 1e-12, 'TolX', 1e-12);
[pFit, ssRes] = fminsearch(cost, p0, opts);
fprintf('Chiusazza fit: a = %.4f, b = %.4f, tau = %.0f ohm*m, SSres = %.2e\n', pFit, ssRes);
rFine = logspace(log10(min(rhoC)*0.9), log10(max(rhoC)*1.05), 300);

figure('Color', 'w', 'Position', [60 60 1300 520]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; box on; grid on;
plot(rhoM, chiM, 'o', 'MarkerFaceColor', cMicio, 'MarkerEdgeColor', 'k', 'MarkerSize', 9, ...
     'DisplayName', 'Constrained model (data)');
yline(freeM, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
      'DisplayName', sprintf('Free model (\\chi^2 = %.3f)', freeM));
set(gca, 'XScale', 'log'); ylim([0.15 0.68]); xlim([1.4e4 1.15e5]);
xlabel('Imposed resistivity, \rho_{input} (\Omega\cdotm)'); ylabel('Reduced \chi^2');
title({'(a) Micio Conti (lava tube)', 'small decrease, no fit (range 0.004)'}, 'FontWeight', 'normal');
legend('Location', 'east');

nexttile; hold on; box on; grid on;
plot(rFine, model(pFit, rFine), '-', 'Color', cChiu, 'LineWidth', 2, 'DisplayName', 'Saturation fit');
plot(rhoC, chiC, 's', 'MarkerFaceColor', cChiu, 'MarkerEdgeColor', 'k', 'MarkerSize', 9, ...
     'DisplayName', 'Constrained model (data)');
yline(freeC, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
      'DisplayName', sprintf('Free model (\\chi^2 = %.3f)', freeC));
set(gca, 'XScale', 'log'); ylim([0.15 0.68]); xlim([3.3e3 1.15e5]);
xlabel('Imposed resistivity, \rho_{input} (\Omega\cdotm)'); ylabel('Reduced \chi^2');
title({'(b) Chiusazza (karst cave)', 'gradual saturation, \tau within tested range'}, 'FontWeight', 'normal');
legend('Location', 'east');

exportgraphics(gcf, 'Figure_chi2_saturation.png', 'Resolution', 300);
