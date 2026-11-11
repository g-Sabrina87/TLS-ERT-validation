%% 7_5_alpha_sensitivity_and_distances.m
% Two checks without new inversions, for Micio Conti and Chiusazza:
%  (A) sensitivity of the result to the definition of the cavity (alpha of the ground
%      truth): AUC of U and H, and gain H-U, for different alpha values of the
%      alpha shape of the full TLS. The rest of the protocol is identical to
%      7_2_null_diagnostics_buffer.m (same evaluation set: held-out rectangle, buffer
%      from the constraint points, within 3 m of the TLS).
%  (B) distance between the constraint region and the held-out sector.
% Output: alpha_sensitivity_<site>.csv

clear; clc;

site = 'Chiusazza';        % 'Micio Conti' or 'Chiusazza'

switch site
    case 'Micio Conti'
        modelFile_U = 'Model_Micio Conti_U.txt';
        modelFile_H = 'Model_Micio Conti_H_rho75000.txt';
        holdoutFile = 'Micio Conti_holdout_region.mat';
        constraintFile = 'Micio Conti_constraint_only.txt';
        tlsFullFile = 'Micio Conti_cropped_ERTarea.txt';
        alphaList = [0.95 1.35 1.45];                 % 1.35 = value used in the paper
        bufferRadius = 1.5;
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93;
    case 'Chiusazza'
        modelFile_U = 'Model_Chiusazza_U.txt';
        modelFile_H = 'Model_Chiusazza_H_rho20000.txt';
        holdoutFile = 'Chiusazza_holdout_region.mat';
        constraintFile = 'Chiusazza_constraint_only.txt';
        tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';
        alphaList = [2.0 2.5 3.5];                     % 2.5 = value used in the paper
        bufferRadius = 1.0;
        polyXY = [ 514191.91  4097666.70;
                   514306.99  4097666.70;
                   514306.99  4097785.75;
                   514191.91  4097785.75 ];
        AOI.Zmin = 79.82;
        AOI.Zmax = 115.71;
end
maxDistTLS = 3;

load(holdoutFile, 'HELDOUT');
tls = readmatrix(tlsFullFile);  tls = tls(:,1:3);  tls = tls(all(isfinite(tls),2),:);
cp  = readmatrix(constraintFile); cp = cp(:,1:3);   cp  = cp(all(isfinite(cp),2),:);

% --- elements inside the held-out rectangle (same selection as 7_2_null_diagnostics_buffer.m)
dU = readmatrix(modelFile_U);  dH = readmatrix(modelFile_H);
sel = @(d) inpolygon(d(:,1),d(:,2),polyXY(:,1),polyXY(:,2)) & d(:,3)>=AOI.Zmin & d(:,3)<=AOI.Zmax & ...
           d(:,1)>=HELDOUT.Xmin & d(:,1)<=HELDOUT.Xmax & ...
           d(:,2)>=HELDOUT.Ymin & d(:,2)<=HELDOUT.Ymax & ...
           d(:,3)>=HELDOUT.Zmin & d(:,3)<=HELDOUT.Zmax;
dU = dU(sel(dU),:);  dH = dH(sel(dH),:);
if size(dU,1) ~= size(dH,1) || max(abs(dU(:,1:3)-dH(:,1:3)),[],'all') > 1e-6
    error('U and H have different elements inside the held-out rectangle.');
end
XYZ = dU(:,1:3);

[~, dC] = knnsearch(cp,  XYZ, 'K', 1);     % distance from the constraint points
[~, dT] = knnsearch(tls, XYZ, 'K', 1);     % distance from the full TLS
keep = (dC > bufferRadius) & (dT <= maxDistTLS);
fprintf('%s: elements inside the held-out rectangle %d, in the evaluation set %d\n', site, size(XYZ,1), sum(keep));

%% (A) sensitivity to the alpha of the ground truth
rows = {};
for a = alphaList
    shp = alphaShape(tls(:,1), tls(:,2), tls(:,3), a);
    truth = inShape(shp, XYZ(keep,1), XYZ(keep,2), XYZ(keep,3));
    nPos = sum(truth); nNeg = sum(~truth);
    aucs = nan(1,2);
    for m = 1:2
        if m == 1, rho = dU(keep,4); else, rho = dH(keep,4); end
        [~,~,g] = unique(rho); cnt = accumarray(g,1); avg = cumsum(cnt) - (cnt-1)/2; rk = avg(g);
        aucs(m) = (sum(rk(truth)) - nPos*(nPos+1)/2) / (nPos*nNeg);
    end
    fprintf('alpha = %.2f | elements inside the true cavity: %d / %d | AUC U = %.4f, H = %.4f, H-U = %+.4f\n', ...
        a, nPos, sum(keep), aucs(1), aucs(2), aucs(2)-aucs(1));
    rows(end+1,:) = {a, sum(keep), nPos, aucs(1), aucs(2), aucs(2)-aucs(1)}; %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'alpha','nEval','nInsideCavity','AUC_U','AUC_H','dAUC_H_minus_U'});
writetable(T, sprintf('alpha_sensitivity_%s.csv', strrep(site,' ','')));

%% (B) distance between the constraint region and the held-out sector
inHOtls = tls(:,1)>=HELDOUT.Xmin & tls(:,1)<=HELDOUT.Xmax & ...
          tls(:,2)>=HELDOUT.Ymin & tls(:,2)<=HELDOUT.Ymax & ...
          tls(:,3)>=HELDOUT.Zmin & tls(:,3)<=HELDOUT.Zmax;
[~, dHO] = knnsearch(cp, tls(inHOtls,:), 'K', 1);
fprintf('\nTLS points inside the held-out sector: %d\n', sum(inHOtls));
fprintf('Distance of these points from the nearest constraint point: min %.2f m, mean %.2f m, median %.2f m, maximum %.2f m\n', ...
    min(dHO), mean(dHO), median(dHO), max(dHO));
fprintf('Evaluation elements: distance from the constraint points min %.2f, mean %.2f, maximum %.2f m\n', ...
    min(dC(keep)), mean(dC(keep)), max(dC(keep)));
