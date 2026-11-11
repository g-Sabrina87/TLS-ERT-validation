%% 7_2_null_diagnostics_buffer.m
% Evaluation of U, H and the null models in the held-out rectangle, excluding the
% elements close to the CONSTRAINT points (inside the held-out rectangle the
% constraint of H was imposed directly). Radius 0 reproduces the evaluation without
% exclusion. Output: null_buffer_<site>_excl<r>_near<R>.csv (one per configuration).

clear; clc;

site = 'Micio Conti';        % 'Micio Conti' or 'Chiusazza'
nearRadii = [Inf 6 3 1.5];   % m: keep only elements within this distance from the full TLS (Inf = no limit)

switch site
    case 'Micio Conti'
        modelFile_U = 'Model_Micio Conti_U.txt';
        modelFile_H = 'Model_Micio Conti_H_rho75000.txt';
        nullPattern = 'model_NULL_%02d.txt';
        holdoutFile = 'Micio Conti_holdout_region.mat';
        constraintFile = 'Micio Conti_constraint_only.txt';
        tlsFullFile = 'Micio Conti_cropped_ERTarea.txt';
        alphaRadius = 1.35;
        classificationThreshold = 48107.1;
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93;
        thresholds = linspace(1000, 60000, 50);
        nNull = 30;
        radii = [0 1.5 3];   % m: buffer radii from the constraint points (the evaluation set of the paper uses 1.5)
    case 'Chiusazza'
        modelFile_U = 'Model_Chiusazza_U.txt';
        modelFile_H = 'Model_Chiusazza_H_rho20000.txt';
        nullPattern = 'Null_%02d.txt';
        holdoutFile = 'Chiusazza_holdout_region.mat';
        constraintFile = 'Chiusazza_constraint_only.txt';
        tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';
        alphaRadius = 2.50;
        classificationThreshold = 5513;
        polyXY = [ 514191.9125  4097666.704;
                   514306.9915  4097666.704;
                   514306.9915  4097785.755;
                   514191.9125  4097785.755 ];
        AOI.Zmin = 79.82597;
        AOI.Zmax = 115.70548;
        thresholds = linspace(1000, 20000, 50);
        nNull = 20;
        radii = [0 1 3];     % m: buffer radii from the constraint points (the evaluation set of the paper uses 1.0)
end

%% ---------- REFERENCE DATA ----------
load(holdoutFile, 'HELDOUT');
tls = readmatrix(tlsFullFile);
shpFull = alphaShape(tls(:,1), tls(:,2), tls(:,3), alphaRadius);
cp = readmatrix(constraintFile); cp = cp(:,1:3);

names = {'U', 'H'}; files = {modelFile_U, modelFile_H};
for k = 1:nNull
    names{end+1} = sprintf('NULL%02d', k); %#ok<SAGROW>
    files{end+1} = sprintf(nullPattern, k); %#ok<SAGROW>
end
nM = numel(names);

%% ---------- READ THE MODELS INSIDE THE HELD-OUT RECTANGLE ----------
D = cell(nM,1);
XYZref = [];
for m = 1:nM
    if ~isfile(files{m}), fprintf('MISSING: %s\n', files{m}); continue; end
    d = readmatrix(files{m}); XYZ = d(:,1:3); rho = d(:,4);
    inTOT = inpolygon(XYZ(:,1), XYZ(:,2), polyXY(:,1), polyXY(:,2)) & ...
            XYZ(:,3) >= AOI.Zmin & XYZ(:,3) <= AOI.Zmax;
    XYZ = XYZ(inTOT,:); rho = rho(inTOT);
    inHO = XYZ(:,1) >= HELDOUT.Xmin & XYZ(:,1) <= HELDOUT.Xmax & ...
           XYZ(:,2) >= HELDOUT.Ymin & XYZ(:,2) <= HELDOUT.Ymax & ...
           XYZ(:,3) >= HELDOUT.Zmin & XYZ(:,3) <= HELDOUT.Zmax;
    XYZ = XYZ(inHO,:); rho = rho(inHO);
    if isempty(XYZref)
        XYZref = XYZ;
        truthRef = inShape(shpFull, XYZ(:,1), XYZ(:,2), XYZ(:,3));
        dC = minDistToSet(XYZ, cp);          % distance of each element from the constraint points
        dT = minDistToSet(XYZ, tls(:,1:3));  % distance of each element from the full TLS
    elseif size(XYZ,1) ~= size(XYZref,1) || max(abs(XYZ(:) - XYZref(:))) > 1e-6
        fprintf('WARNING: %s has different elements from the first model\n', names{m});
    end
    D{m} = rho;
end

fprintf('\nElements inside the held-out rectangle: %d\n', size(XYZref,1));
for r = radii
    fprintf('  radius %.1f m: elements kept %d (excluded %d)\n', r, sum(dC > r), sum(dC <= r));
end

%% ---------- METRICS FOR EACH RADIUS ----------
for r = radii
  for R = nearRadii
    keep = (dC > r) & (dT <= R);
    truth = truthRef(keep);
    res = nan(nM, 4);   % IoU, AUC (ranks), nAboveThr, grid AUC
    for m = 1:nM
        if isempty(D{m}), continue; end
        rho = D{m}(keep);
        pred = rho >= classificationThreshold;
        TP = sum(pred & truth); FP = sum(pred & ~truth); FN = sum(~pred & truth);
        nPos = sum(truth); nNeg = sum(~truth);
        if nPos > 0 && nNeg > 0
            % rank AUC (exact, no threshold grid)
            [~,~,g] = unique(rho); cnt = accumarray(g,1); avg = cumsum(cnt) - (cnt-1)/2; ranks = avg(g);
            auc = (sum(ranks(truth)) - nPos*(nPos+1)/2) / (nPos*nNeg);
            % AUC on a threshold grid (for comparison)
            tpr = arrayfun(@(t) sum(rho >= t & truth)/nPos, thresholds);
            fpr = arrayfun(@(t) sum(rho >= t & ~truth)/nNeg, thresholds);
            [fs, ix] = sort(fpr); aucGrid = abs(trapz(fs, tpr(ix)));
        else
            auc = NaN; aucGrid = NaN;
        end
        res(m,:) = [TP/max(TP+FP+FN,1), auc, sum(pred), aucGrid];
    end
    % values of each model for this configuration
    T2 = table(string(names(:)), res(:,2), res(:,1), res(:,3), res(:,4), ...
        'VariableNames', {'model','AUC_rank','IoU','nAboveThr','AUC_grid'});
    writetable(T2, sprintf('null_buffer_%s_excl%g_near%g.csv', regexprep(site,'\s+','_'), r, R));
    iouU = res(1,1); iouH = res(2,1); aucU = res(1,2); aucH = res(2,2);
    iouN = res(3:end,1); aucN = res(3:end,2);
    iouN = iouN(~isnan(iouN)); aucN = aucN(~isnan(aucN));
    fprintf('\n===== %s, buffer %.1f m, only within %.1f m of the TLS: %d elements, %d inside the true cavity =====\n', ...
        site, r, R, sum(keep), sum(truth));
    fprintf('AUC: U=%.4f H=%.4f | nulls: mean=%.4f (min %.4f, max %.4f)\n', aucU, aucH, mean(aucN), min(aucN), max(aucN));
    fprintf('     nulls with AUC > U: %d of %d ; nulls with AUC >= H: %d of %d\n', sum(aucN>aucU), numel(aucN), sum(aucN>=aucH), numel(aucN));
    fprintf('IoU: U=%.4f H=%.4f | nulls: mean=%.4f (min %.4f, max %.4f)\n', iouU, iouH, mean(iouN), min(iouN), max(iouN));
    fprintf('     nulls with IoU >= H: %d of %d\n', sum(iouN>=iouH), numel(iouN));
    fprintf('(check) threshold-grid AUC: U=%.4f H=%.4f ; ranks: U=%.4f H=%.4f\n', res(1,4), res(2,4), aucU, aucH);
    fprintf('AUC gain over U: H=%+.4f ; nulls (mean)=%+.4f\n', aucH-aucU, mean(aucN)-aucU);
  end
end

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
