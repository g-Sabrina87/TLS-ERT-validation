%% 7_4_2_iou_threshold_sweep_Chiusazza.m
% As 7_4_1_iou_threshold_sweep_MicioConti.m, for Chiusazza.
% IoU of U, H and the nulls over a RANGE of thresholds, on the same evaluation set
% as the paper (11,758 elements, 2,200 in the cavity): held-out rectangle, farther
% than 1.0 m from the constraint points, within 3 m of the full TLS. Mean IoU,
% maximum IoU and Average Precision (threshold-free) are also computed, and the same
% measures for the nulls with little overlap with the TLS (< farLimit % of the points
% within 3 m).
% Output: iou_sweep_Chiusazza.csv, iou_sweep_summary_Chiusazza.csv,
%         Figure_iou_sweep_Chiusazza.png

clear; clc; close all;

%% ---------- SETTINGS ----------
modelFile_U = 'Model_Chiusazza_U.txt';
modelFile_H = 'Model_Chiusazza_H_rho20000.txt';
nullPattern = 'Null_%02d.txt';
holdoutFile = 'Chiusazza_holdout_region.mat';
constraintFile = 'Chiusazza_constraint_only.txt';
tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';
alphaRadius = 2.50;

polyXY = [ 514191.9125  4097666.704;
           514306.9915  4097666.704;
           514306.9915  4097785.755;
           514191.9125  4097785.755 ];
AOI.Zmin = 79.82597;
AOI.Zmax = 115.70548;

bufferRadius = 1.0;      % exclusion around the constraint points (m)
maxDistTLS   = 3;        % keep only elements within this distance from the TLS (m)
refThreshold = 5513;     % threshold of the paper for Chiusazza

thresholds = logspace(log10(1500), log10(40000), 80);   % threshold range
nNull = 20;

% Overlap of each null with the full TLS (% of points within 3 m), read from the output
% of 4_4_null_position_check.m (site = 'Chiusazza'), column pctNearCavity; row NullID = 0
% is the real constraint and is dropped. Used only for the subset of "far" nulls and
% for the rank correlation with the average precision.
positionFile = 'null_position_Chiusazza.csv';
if ~isfile(positionFile)
    error('File not found: %s (output of 4_4_null_position_check.m). Current folder: %s', positionFile, pwd);
end
Ppos = readtable(positionFile);
Ppos = sortrows(Ppos(Ppos.NullID > 0, :), 'NullID');
if height(Ppos) ~= nNull
    error('%s has %d nulls, expected %d.', positionFile, height(Ppos), nNull);
end
overlap3m = Ppos.pctNearCavity(:)';     % % of null points within 3 m of the full TLS, NULL01 ... NULL20
farLimit = 15;           % "far" nulls: overlap < farLimit %

%% ---------- REFERENCE DATA ----------
load(holdoutFile, 'HELDOUT');
tls = readmatrix(tlsFullFile);
shpFull = alphaShape(tls(:,1), tls(:,2), tls(:,3), alphaRadius);
cp = readmatrix(constraintFile); cp = cp(:,1:3);

names = {'U','H'}; files = {modelFile_U, modelFile_H};
for k = 1:nNull
    names{end+1} = sprintf('NULL%02d',k); %#ok<SAGROW>
    files{end+1} = sprintf(nullPattern,k); %#ok<SAGROW>
end
nM = numel(names);

%% ---------- READ THE MODELS ----------
D = cell(nM,1); XYZref = [];
for m = 1:nM
    if ~isfile(files{m}), fprintf('MISSING: %s\n', files{m}); continue; end
    d = readmatrix(files{m}); XYZ = d(:,1:3); rho = d(:,4);
    inTOT = inpolygon(XYZ(:,1),XYZ(:,2),polyXY(:,1),polyXY(:,2)) & ...
            XYZ(:,3) >= AOI.Zmin & XYZ(:,3) <= AOI.Zmax;
    XYZ = XYZ(inTOT,:); rho = rho(inTOT);
    inHO = XYZ(:,1) >= HELDOUT.Xmin & XYZ(:,1) <= HELDOUT.Xmax & ...
           XYZ(:,2) >= HELDOUT.Ymin & XYZ(:,2) <= HELDOUT.Ymax & ...
           XYZ(:,3) >= HELDOUT.Zmin & XYZ(:,3) <= HELDOUT.Zmax;
    XYZ = XYZ(inHO,:); rho = rho(inHO);
    if isempty(XYZref)
        XYZref = XYZ;
        truthRef = inShape(shpFull, XYZ(:,1), XYZ(:,2), XYZ(:,3));
        dC = minDistToSet(XYZ, cp);
        dT = minDistToSet(XYZ, tls(:,1:3));
    elseif size(XYZ,1) ~= size(XYZref,1) || max(abs(XYZ(:)-XYZref(:))) > 1e-6
        fprintf('WARNING: %s has different elements from the first model\n', names{m});
    end
    D{m} = rho;
end

fprintf('Elements inside the held-out rectangle %d (expected 87870); with dC>%.1f m: %d (expected 87752); with dT<=%.1f m: %d\n', ...
    size(XYZref,1), bufferRadius, sum(dC>bufferRadius), maxDistTLS, sum(dT<=maxDistTLS));
fprintf('Constraint points %d, TLS points %d, U file = %s\n', size(cp,1), size(tls,1), files{1});
keep  = (dC > bufferRadius) & (dT <= maxDistTLS);
truth = truthRef(keep);
nPos = sum(truth); nNeg = sum(~truth);
fprintf('Evaluation set: %d elements, %d inside the true cavity (expected: 11758 and 2200)\n', ...
    sum(keep), nPos);

%% ---------- IoU vs THRESHOLD + AP ----------
nT = numel(thresholds);
IoU  = nan(nT, nM);
nAbv = nan(nT, nM);
AP   = nan(1, nM);
for m = 1:nM
    if isempty(D{m}), continue; end
    rho = D{m}(keep);
    for i = 1:nT
        pred = rho >= thresholds(i);
        TP = sum(pred & truth); FP = sum(pred & ~truth); FN = sum(~pred & truth);
        IoU(i,m)  = TP / max(TP+FP+FN, 1);
        nAbv(i,m) = sum(pred);
    end
    % Average Precision without threshold (vectorized, handles ties)
    [rs, ord] = sort(rho, 'descend'); tt = truth(ord);
    tp = cumsum(tt); kk = (1:numel(rs))';
    lastTie = [rs(1:end-1) ~= rs(2:end); true];
    P = tp(lastTie) ./ kk(lastTie); Rc = tp(lastTie) / nPos;
    AP(m) = sum(diff([0; Rc]) .* P);
end

%% ---------- SUMMARY ----------
logT = log(thresholds(:));
span = logT(end) - logT(1);
meanIoU = trapz(logT, IoU) / span;          % mean IoU over the range (log scale)
maxIoU  = max(IoU, [], 1);
[~, iRef] = min(abs(thresholds - refThreshold));

iU = 1; iH = 2; iN = 3:nM;
dIoU = IoU(:,iH) - IoU(:,iU);               % H - U at each threshold
fprintf('\n===== IoU H - U over the threshold range %.0f - %.0f Ohm*m =====\n', thresholds(1), thresholds(end));
fprintf('Thresholds with H < U: %d of %d ; H = U: %d ; H > U: %d\n', ...
    sum(dIoU<-1e-9), nT, sum(abs(dIoU)<=1e-9), sum(dIoU>1e-9));
fprintf('Difference H - U: mean %+.3f, min %+.3f, max %+.3f\n', mean(dIoU,'omitnan'), min(dIoU), max(dIoU));
fprintf('At the threshold of the paper (~%.0f): IoU U=%.3f H=%.3f (above threshold: U=%d, H=%d)\n', ...
    thresholds(iRef), IoU(iRef,iU), IoU(iRef,iH), nAbv(iRef,iU), nAbv(iRef,iH));
fprintf('\nMean IoU over the range: U=%.3f H=%.3f ; nulls: mean %.3f (min %.3f, max %.3f)\n', ...
    meanIoU(iU), meanIoU(iH), mean(meanIoU(iN),'omitnan'), min(meanIoU(iN)), max(meanIoU(iN)));
fprintf('Maximum IoU over the range: U=%.3f H=%.3f ; nulls: mean %.3f (min %.3f, max %.3f)\n', ...
    maxIoU(iU), maxIoU(iH), mean(maxIoU(iN),'omitnan'), min(maxIoU(iN)), max(maxIoU(iN)));
fprintf('Average Precision (no threshold): U=%.3f H=%.3f ; nulls: mean %.3f (min %.3f, max %.3f) ; positives %d/%d = %.3f (chance level)\n', ...
    AP(iU), AP(iH), mean(AP(iN),'omitnan'), min(AP(iN)), max(AP(iN)), nPos, numel(truth), nPos/numel(truth));
fprintf('Nulls with mean IoU >= H: %d of %d ; nulls with AP >= H: %d of %d\n', ...
    sum(meanIoU(iN) >= meanIoU(iH)), numel(iN), sum(AP(iN) >= AP(iH)), numel(iN));

% --- nulls "far" from the cavity (overlap < farLimit %) ---
farIdx = iN(overlap3m < farLimit);
fprintf('\nFar nulls (overlap < %g%%): %d (NULL%s)\n', farLimit, numel(farIdx), strjoin(arrayfun(@(k) sprintf('%02d',k-2), farIdx, 'UniformOutput', false), ','));
fprintf('  Mean IoU over the range: far nulls mean %.4f (min %.4f, max %.4f) ; H=%.4f ; U=%.4f\n', ...
    mean(meanIoU(farIdx)), min(meanIoU(farIdx)), max(meanIoU(farIdx)), meanIoU(iH), meanIoU(iU));
fprintf('  AP: far nulls mean %.4f (min %.4f, max %.4f) ; H=%.4f ; U=%.4f\n', ...
    mean(AP(farIdx)), min(AP(farIdx)), max(AP(farIdx)), AP(iH), AP(iU));
fprintf('  far nulls with mean IoU >= H: %d of %d ; with AP >= H: %d of %d\n', ...
    sum(meanIoU(farIdx) >= meanIoU(iH)), numel(farIdx), sum(AP(farIdx) >= AP(iH)), numel(farIdx));
apN = AP(iN); [rhoOv, pOv] = corr(overlap3m(:), apN(:), 'Type', 'Spearman');
fprintf('  Spearman overlap-AP over the 20 nulls: rho=%.2f (p=%.2f, descriptive)\n', rhoOv, pOv);
iouN = meanIoU(iN);
[rhoOvI, pOvI] = corr(overlap3m(:), iouN(:), 'Type', 'Spearman');
fprintf('  Spearman overlap-mean IoU over the 20 nulls: rho=%.2f (p=%.2f, descriptive)\n', rhoOvI, pOvI);
fprintf('  Nulls with AP >= H: NULL%s\n', strjoin(arrayfun(@(k) sprintf('%02d',k-2), iN(AP(iN) >= AP(iH)), 'UniformOutput', false), ',NULL'));

%% ---------- OUTPUT FILES ----------
T = table(thresholds(:), IoU(:,iU), IoU(:,iH), mean(IoU(:,iN),2,'omitnan'), ...
    min(IoU(:,iN),[],2), max(IoU(:,iN),[],2), nAbv(:,iU), nAbv(:,iH), ...
    'VariableNames', {'threshold','IoU_U','IoU_H','IoU_nulls_mean','IoU_nulls_min','IoU_nulls_max','nAbove_U','nAbove_H'});
writetable(T, 'iou_sweep_Chiusazza.csv');

S = table(string(names(:)), meanIoU(:), maxIoU(:), AP(:), ...
    'VariableNames', {'model','IoU_mean_over_thresholds','IoU_max_over_thresholds','AP'});
writetable(S, 'iou_sweep_summary_Chiusazza.csv');

%% ---------- FIGURE ----------
figure('Color','w','Position',[100 100 900 420]);
hold on;
nullMin = min(IoU(:,iN),[],2); nullMax = max(IoU(:,iN),[],2);
fill([thresholds(:); flipud(thresholds(:))], [nullMin; flipud(nullMax)], [0.8 0.8 0.8], ...
    'EdgeColor','none','FaceAlpha',0.6);
plot(thresholds, mean(IoU(:,iN),2,'omitnan'), '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
plot(thresholds, IoU(:,iU), '-',  'Color',[0.00 0.45 0.74], 'LineWidth',2);
plot(thresholds, IoU(:,iH), '-',  'Color',[0.85 0.33 0.10], 'LineWidth',2);
xline(refThreshold, 'k--');
set(gca,'XScale','log'); grid on; box on;
xlabel('Classification threshold (\Omega\cdotm)'); ylabel('IoU');
legend({'Nulls (min-max)','Nulls (mean)','U','H','Paper threshold'}, 'Location','best');
title(sprintf('Chiusazza: IoU vs threshold (%d elements, %d in cavity)', sum(keep), nPos));
exportgraphics(gcf, 'Figure_iou_sweep_Chiusazza.png', 'Resolution', 300);

%% ---------- LOCAL FUNCTION ----------
function dmin = minDistToSet(P, Q)
% minimum distance of each point of P from the set Q (in blocks)
    o = mean(Q, 1); P = P - o; Q = Q - o;
    q2 = sum(Q.^2, 2)';
    dmin = zeros(size(P,1), 1);
    chunk = 200;
    for i = 1:chunk:size(P,1)
        ii = i:min(i+chunk-1, size(P,1));
        Pc = P(ii,:);
        D2 = sum(Pc.^2, 2) + q2 - 2*(Pc*Q');
        dmin(ii) = sqrt(max(min(D2, [], 2), 0));
    end
end
