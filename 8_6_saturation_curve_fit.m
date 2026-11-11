%% 8_6_saturation_curve_fit.m
% Fit of an exponential saturation curve to the reduced chi-square as a function of the
% resistivity imposed at the constraint:
%
%       chi2(rho) = a + b * exp(-rho / tau)
%
%  a   : asymptotic value (plateau) for very strong constraints (rho -> infinity)
%  b   : amplitude of the effect
%  tau : characteristic saturation scale
%
% The fit uses ONLY the constrained points (the free model is not the limit case
% rho_input -> 0, but a qualitatively different model); the free value is shown as a
% reference line. Uses fminsearch (base MATLAB, no toolbox required).

clear; clc; close all;

%% ---------- 1) DATA ----------
siteName = 'Chiusazza';

% Constrained points: imposed resistivity and corresponding reduced chi-square
rho_data  = [4000, 10000, 20000, 50000, 100000];
chi2_data = [0.3666, 0.347, 0.3391, 0.3326, 0.3239];

% Free (unconstrained) model, reference only, not included in the fit
chi2_free = 0.4188;

%% ---------- 2) FIT ----------
model = @(p, rho) p(1) + p(2) * exp(-rho ./ p(3));
costFun = @(p) sum((model(p, rho_data) - chi2_data).^2);

a0   = min(chi2_data);
b0   = chi2_data(1) - a0;
tau0 = mean(rho_data) / 3;
p0   = [a0, b0, tau0];

options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'TolFun', 1e-10, 'TolX', 1e-10);
[pFit, resnorm] = fminsearch(costFun, p0, options);

a_fit   = pFit(1);
b_fit   = pFit(2);
tau_fit = pFit(3);

fprintf('=== Saturation curve fit - %s ===\n', siteName);
fprintf('Model: chi2(rho) = a + b*exp(-rho/tau)\n\n');
fprintf('a   (asymptotic plateau)  = %.4f\n', a_fit);
fprintf('b   (effect amplitude)    = %.4f\n', b_fit);
fprintf('tau (saturation scale)    = %.1f ohm*m\n', tau_fit);
fprintf('Sum of squared residuals  = %.6f\n\n', resnorm);

chi2_predicted = model(pFit, rho_data);
fprintf('Observed vs predicted:\n');
for i = 1:numel(rho_data)
    fprintf('  rho=%6d -> observed=%.4f, predicted=%.4f, residual=%+.4f\n', ...
        rho_data(i), chi2_data(i), chi2_predicted(i), chi2_data(i)-chi2_predicted(i));
end

ssRes = sum((chi2_data - chi2_predicted).^2);
ssTot = sum((chi2_data - mean(chi2_data)).^2);
Rsq = 1 - ssRes/ssTot;
fprintf('\nR^2 (N=%d points): %.4f\n', numel(rho_data), Rsq);

%% ---------- 3) FIGURE ----------
rhoFine = linspace(0, max(rho_data)*1.3, 500);
chi2Fine = model(pFit, rhoFine);

figure;
hold on;
plot(rhoFine, chi2Fine, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Exponential fit');
plot(rho_data, chi2_data, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, ...
    'DisplayName', 'Observed (constrained)');
yline(chi2_free, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Free (%.4f)', chi2_free));
yline(a_fit, 'b:', 'LineWidth', 1.2, 'Label', sprintf('Estimated plateau (%.4f)', a_fit));
xlabel('Resistivity imposed at the constraint, rho\_input (\Omega\cdotm)');
ylabel('Reduced chi-square');
title(sprintf('Chi-square saturation curve - %s', siteName));
legend('Location', 'best');
grid on;
set(gca, 'XScale', 'linear');
