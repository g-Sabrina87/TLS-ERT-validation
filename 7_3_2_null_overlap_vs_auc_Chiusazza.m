%% 7_3_2_null_overlap_vs_auc_Chiusazza.m
% Does the AUC of the null constraints depend on how much of each null overlaps (a) the
% held-out sector and (b) the cavity (full TLS cloud)? At Chiusazza all 20 nulls touch the
% TLS (minimum distance < 0.1 m), so the distance used at Micio Conti (7_3) does not vary;
% the quantities that vary are the shares of null points inside the held-out sector and
% within 3 m of the TLS. Pearson and Spearman (1904) correlations with the AUC of the
% nulls on the main evaluation set (exclusion 1 m from the constraint points, elements
% within 3 m of the full TLS cloud). Descriptive only: the nulls are translated copies of
% one shape and are not independent replicates. No inversions.
%
% Inputs (Chiusazza, 20 nulls):
%   null_position_Chiusazza.csv                 from 4_4_null_position_check.m (site = 'Chiusazza')
%                                               columns: NullID, tx, ty, tz, pctInHeldout, pctNearCavity,
%                                               minDistTLS, medDistElec, nElec10 (row NullID = 0 is the real constraint)
%   null_buffer_Chiusazza_excl1_near3.csv       from 7_2_null_diagnostics_buffer.m (site = 'Chiusazza')
%                                               columns: model, AUC_rank, IoU, nAboveThr, AUC_grid
% Output: printed summary (and null_overlap_vs_auc_Chiusazza.csv).

clear; clc;

positionFile = 'null_position_Chiusazza.csv';
aucFile      = 'null_buffer_Chiusazza_excl1_near3.csv';

for f = {positionFile, aucFile}
    if ~isfile(f{1})
        error(['File not found: %s\nCurrent folder: %s\n' ...
               'Put the two input files in the current folder or edit the names above.'], f{1}, pwd);
    end
end

P = readtable(positionFile);
P = P(P.NullID > 0, :);                 % drop row 0 (the real constraint, zero translation)
T = readtable(aucFile);

aucU = T.AUC_rank(strcmp(T.model, 'U'));
aucH = T.AUC_rank(strcmp(T.model, 'H'));

aucN = nan(height(P), 1);
for k = 1:height(P)
    aucN(k) = T.AUC_rank(strcmp(T.model, sprintf('NULL%02d', P.NullID(k))));
end
if any(isnan(aucN))
    error('Some nulls of %s are missing from %s.', positionFile, aucFile);
end

dAUC = aucN - aucU;
inHO = P.pctInHeldout;      % % of null points inside the held-out sector
near = P.pctNearCavity;     % % of null points within 3 m of the full TLS cloud

[rP_ho, pP_ho] = corr(inHO, aucN, 'Type', 'Pearson');
[rS_ho, pS_ho] = corr(inHO, aucN, 'Type', 'Spearman');
[rP_nc, pP_nc] = corr(near, aucN, 'Type', 'Pearson');
[rS_nc, pS_nc] = corr(near, aucN, 'Type', 'Spearman');

fprintf('\nChiusazza, %d nulls\n', numel(aucN));
fprintf('  AUC: U = %.4f, H = %.4f\n', aucU, aucH);
fprintf('  Nulls: mean %.4f, range %.4f-%.4f; with AUC >= H: %d of %d\n', ...
    mean(aucN), min(aucN), max(aucN), sum(aucN >= aucH), numel(aucN));
fprintf('  Mean AUC difference from U: nulls %+.4f, H %+.4f\n', mean(dAUC), aucH - aucU);
fprintf('  Share of null points inside the held-out sector: mean %.1f%% (%.1f-%.1f%%); >10%%: %d of %d\n', ...
    mean(inHO), min(inHO), max(inHO), sum(inHO > 10), numel(inHO));
fprintf('  Share of null points within 3 m of the TLS: median %.1f%% (%.1f-%.1f%%)\n', ...
    median(near), min(near), max(near));
fprintf('  AUC vs share inside the held-out sector: Pearson r = %+.3f (p = %.2f), Spearman rho = %+.2f (p = %.2f)\n', ...
    rP_ho, pP_ho, rS_ho, pS_ho);
fprintf('  AUC vs share within 3 m of the TLS:      Pearson r = %+.3f (p = %.2f), Spearman rho = %+.2f (p = %.2f)\n', ...
    rP_nc, pP_nc, rS_nc, pS_nc);
fprintf('  (descriptive: the nulls are not independent replicates)\n');

R = table(P.NullID, inHO, near, aucN, dAUC, ...
    'VariableNames', {'NullID','pctInHeldout','pctNearCavity','AUC','dAUC_from_U'});
writetable(R, 'null_overlap_vs_auc_Chiusazza.csv');
