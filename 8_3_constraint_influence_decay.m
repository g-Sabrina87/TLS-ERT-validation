%% 8_3_constraint_influence_decay.m
% Distance over which the constraint changes the inverted model. On already inverted
% models, computes the difference of log10(resistivity) between H and U as a function
% of the distance from the constraint cloud, and compares it with that of the null
% models (null - U as a function of the distance from the null's own translated cloud).
% The last distance class is the far-field level of the difference between two
% different inversions.
% Output: influence_decay_<site>.csv and influence_decay_<site>.png

clear; clc;

site = 'Chiusazza';        % 'Micio Conti' or 'Chiusazza'
useNulls = true;           % false to skip the null models
edges = [0 1 2 3 5 8 12 20 Inf];     % distance classes from the constraint (m)
nullCloudDir = ['NULL_TLS_constraints_' strrep(site,' ','')];   % folder of the translated clouds written by 4_1 ('' = current folder)

switch site
    case 'Chiusazza'
        modelFile_U = 'Model_Chiusazza_U.txt';
        modelFile_H = 'Model_Chiusazza_H_rho20000.txt';
        nullModelPattern = 'Null_%02d.txt';
        nullCloudPattern = 'TLS_NULL_%02d.txt';          % constraint clouds of the nulls
        nNull = 20; excludeNulls = [7 10 20];            % nulls with a different iteration sequence (as in Table 5)
        constraintFile = 'Chiusazza_constraint_only.txt';
        polyXY = [ 514191.91  4097666.70;
                   514306.99  4097666.70;
                   514306.99  4097785.75;
                   514191.91  4097785.75 ];
        AOI.Zmin = 79.82;
        AOI.Zmax = 115.71;
    case 'Micio Conti'
        modelFile_U = 'Model_Micio Conti_U.txt';
        modelFile_H = 'Model_Micio Conti_H_rho75000.txt';
        nullModelPattern = 'Model_Null_%02d.txt';
        nullCloudPattern = 'TLS_NULL_%02d.txt';
        nNull = 30; excludeNulls = [];
        constraintFile = 'Micio Conti_constraint_only.txt';
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93;
end

% ---------- FILE CHECK ----------
req = {modelFile_U, modelFile_H, constraintFile};
miss = req(~cellfun(@isfile, req));
if ~isempty(miss), error('Files not found in the current folder (%s):\n  %s', pwd, strjoin(miss, sprintf('\n  '))); end

%% ---------- U ----------
dU = readmatrix(modelFile_U);
XYZU = dU(:,1:3); rhoU = dU(:,4);
inA = inpolygon(XYZU(:,1),XYZU(:,2),polyXY(:,1),polyXY(:,2)) & XYZU(:,3)>=AOI.Zmin & XYZU(:,3)<=AOI.Zmax;
if strcmp(site,'Micio Conti')
    inA = true(size(XYZU,1),1);   % whole mesh at Micio Conti (the nulls lie outside the area of interest)
end
P = XYZU(inA,:); logU = log10(rhoU(inA));
fprintf('%s: %d elements in the evaluated area\n', site, size(P,1));
nB = numel(edges)-1;
binOf = @(dist) discretize(dist, edges);

%% ---------- H ----------
cp = readmatrix(constraintFile); cp = cp(all(isfinite(cp(:,1:3)),2),1:3);
dLogH = deltaLog(modelFile_H, XYZU, inA, logU);
[~, dC] = knnsearch(cp, P, 'K', 1);
bH = binOf(dC);
H_med = nan(nB,1); H_mean = nan(nB,1); H_p90 = nan(nB,1); nEl = zeros(nB,1);
for b = 1:nB
    s = (bH == b); nEl(b) = nnz(s);
    if nEl(b) > 0
        H_med(b) = median(abs(dLogH(s))); H_mean(b) = mean(dLogH(s)); H_p90(b) = prctile(abs(dLogH(s)),90);
    end
end

%% ---------- NULL MODELS ----------
N_med = nan(nB, nNull);
if useNulls
    for k = 1:nNull
        if ismember(k, excludeNulls), continue; end
        fM = sprintf(nullModelPattern,k); fC = fullfile(nullCloudDir, sprintf(nullCloudPattern,k));
        if ~isfile(fM) || ~isfile(fC), fprintf('Null %02d: missing file (%s / %s), skipped\n',k,fM,fC); continue; end
        dl = deltaLog(fM, XYZU, inA, logU);
        Q = readmatrix(fC); Q = Q(all(isfinite(Q(:,1:3)),2),1:3);
        [~, dQ] = knnsearch(Q, P, 'K', 1);
        bq = binOf(dQ);
        for b = 1:nB
            s = (bq == b); if nnz(s) > 0, N_med(b,k) = median(abs(dl(s))); end
        end
    end
end
% nFacMean: mean of the per-null factors
nMean = mean(N_med,2,'omitnan'); nFacMean = mean(10.^N_med,2,'omitnan'); nMin = min(N_med,[],2,'omitnan'); nMax = max(N_med,[],2,'omitnan');

%% ---------- OUTPUT ----------
lab = arrayfun(@(i) sprintf('%g-%g', edges(i), edges(i+1)), 1:nB, 'UniformOutput', false)';
T = table(lab, nEl, H_med, 10.^H_med, H_mean, H_p90, nMean, 10.^nMean, nFacMean, nMin, nMax, ...
    'VariableNames', {'dist_m','nElem','H_median_absDlog','H_factor','H_mean_Dlog','H_p90_absDlog','nulls_mean_of_medians','nulls_factor','nulls_meanOfFactors','nulls_min','nulls_max'});
disp(T); writetable(T, sprintf('influence_decay_%s.csv', strrep(site,' ','')));
fprintf('\nDelta = log10(rho_H) - log10(rho_U) (1 unit = factor 10). Far field (last class): H %.3f ; nulls %.3f.\n', H_med(end), nMean(end));

figure; hold on; xm = 1:nB;
if useNulls, fill([xm fliplr(xm)], [nMin' fliplr(nMax')], [0.85 0.85 0.85], 'EdgeColor','none'); plot(xm, nMean, 'k-o'); end
plot(xm, H_med, 'r-s', 'LineWidth', 1.5);
set(gca,'XTick',xm,'XTickLabel',lab); xlabel('Distance from the constraint cloud (m)'); ylabel('Median of |\Delta log_{10}\rho| relative to U');
if useNulls, legend({'Nulls: min-max','Nulls: mean','H'},'Location','northeast'); else, legend({'H'}); end
title(site); grid on;
saveas(gcf, sprintf('influence_decay_%s.png', strrep(site,' ','')));

%% ---------- LOCAL FUNCTION ----------
function dl = deltaLog(modelFile, XYZU, inA, logU)
    d = readmatrix(modelFile); X = d(:,1:3); r = d(:,4);
    if size(X,1) == size(XYZU,1) && max(abs(X(:)-XYZU(:))) < 1e-6
        rr = r(inA);
    else
        [ix, dd] = knnsearch(X, XYZU(inA,:), 'K', 1);
        if max(dd) > 0.01, error('Mesh different from that of U in %s (max distance %.3f m)', modelFile, max(dd)); end
        rr = r(ix);
    end
    dl = log10(rr) - logU;
end
