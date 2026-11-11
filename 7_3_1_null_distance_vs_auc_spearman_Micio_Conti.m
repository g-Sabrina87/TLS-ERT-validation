%% 7_3_null_distance_vs_auc_spearman.m
% Does the AUC of the null constraints depend on their distance from the cavity?
% Spearman rank correlation (Spearman, 1904) between the minimum distance of each null
% from the full TLS cloud and the absolute AUC difference of the null from U, on the
% main evaluation set (exclusion 1.5 m from the constraint points, elements within 3 m
% of the full TLS cloud). Descriptive only: the nulls are translated copies of one shape
% and are not independent replicates. No inversions.
%
% Inputs (Micio Conti, second series of 30 nulls):
%   null_buffer_Micio_Conti_excl1.5_near3.csv   from 7_2_null_diagnostics_buffer.m
%                                               columns: model, AUC_rank, IoU, nAboveThr, AUC_grid
%   NULL_summary_Micio_Conti.csv                NULL_summary_startNN.csv written by 4_1_make_null_constraint_from_TLS.m
%                                               (renamed for Micio Conti); columns: NullID, dx_m, dy_m, dz_m, MinDistToTLS_m, ...
% Output: printed summary (and null_distance_vs_auc_Micio_Conti.csv).

clear; clc;

aucFile     = 'null_buffer_Micio_Conti_excl1.5_near3.csv';
summaryFile = 'NULL_summary_Micio Conti.csv';
nearestK    = 6;      % number of nearest nulls reported separately

T = readtable(aucFile);
S = readtable(summaryFile);

aucU = T.AUC_rank(strcmp(T.model, 'U'));
aucH = T.AUC_rank(strcmp(T.model, 'H'));

% AUC of each null, matched to the NullID of the summary file (NULL01 <-> NullID 1)
aucN = nan(height(S), 1);
for k = 1:height(S)
    aucN(k) = T.AUC_rank(strcmp(T.model, sprintf('NULL%02d', S.NullID(k))));
end
if any(isnan(aucN))
    error('Some nulls of the summary file are missing from %s.', aucFile);
end

dAUC = aucN - aucU;
dist = S.MinDistToTLS_m;

[rho, p] = corr(dist, abs(dAUC), 'Type', 'Spearman');   % ties handled with mean ranks

[~, order] = sort(dist);
near = order(1:nearestK);

fprintf('\nMicio Conti, %d nulls\n', numel(aucN));
fprintf('  AUC: U = %.4f, H = %.4f\n', aucU, aucH);
fprintf('  Nulls: range %.4f-%.4f; with AUC > U: %d of %d; with AUC >= H: %d of %d\n', ...
    min(aucN), max(aucN), sum(aucN > aucU), numel(aucN), sum(aucN >= aucH), numel(aucN));
fprintf('  Mean AUC difference from U: nulls %+.4f, H %+.4f\n', mean(dAUC), aucH - aucU);
fprintf('  Minimum distance from the TLS: %.1f-%.1f m (mean %.1f m)\n', min(dist), max(dist), mean(dist));
fprintf('  Spearman rho (distance vs |dAUC|) = %.2f (p = %.2f, descriptive)\n', rho, p);
fprintf('  %d nearest nulls (%.1f-%.1f m): maximum |dAUC| from U = %.4f\n', ...
    nearestK, min(dist(near)), max(dist(near)), max(abs(dAUC(near))));

R = table(S.NullID, dist, aucN, dAUC, 'VariableNames', {'NullID','MinDistToTLS_m','AUC','dAUC_from_U'});
writetable(R, 'null_distance_vs_auc_Micio_Conti.csv');
