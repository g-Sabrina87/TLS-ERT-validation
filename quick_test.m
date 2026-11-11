%% quick_test.m
% Quick test of the repository: Micio Conti, main evaluation set.
%
% What it does (about 1-2 minutes, no ERTLab64 needed):
%   1. splits the TLS cloud into constraint region and held-out sector;
%   2. builds the main evaluation set (mesh elements inside the held-out
%      rectangle, farther than 1.5 m from the constraint points, within 3 m
%      of the full TLS cloud) and the ground truth (alpha shape of the full
%      TLS cloud, alpha = 1.35 m);
%   3. computes the rank AUC of the inverted models U, H and C on it;
%   4. runs the paired spatial block bootstrap (blocks 3 x 3 x 2 m) for
%      H - U and C - U;
%   5. compares the results with the values of the paper and prints PASS/FAIL.
%
% Requirements: MATLAB R2025b (earlier releases from R2018b should work),
% Statistics and Machine Learning Toolbox (prctile in older releases).
%
% Input files (searched anywhere under the folder of this script):
%   - Micio Conti_cropped_ERTarea.txt   full TLS cloud (65,742 points)
%   - Model_Micio Conti_U.txt, Model_Micio Conti_H_rho75000.txt,
%     Model_Micio Conti_C_rho75000.txt  inverted models (X Y Z resistivity ...)
% Run:  >> quick_test

clear; clc;
tStart = tic;
root = fileparts(mfilename('fullpath'));

%% 1) Settings: values of the paper (Micio Conti)
alphaRadius  = 1.35;      % m, alpha-shape radius of the ground truth
bufferRadius = 1.5;       % m, exclusion distance from the constraint points
maxDistTLS   = 3;         % m, maximum distance from the full TLS cloud
blockSize    = [3 3 2];   % m, bootstrap block (x, y, z)
nBootstrap   = 2000;      % replicates
rng(1);

% ERT quadrilateral and depth range of the model
polyERT = [ 510340.5  4157110;
            510311.1  4157106;
            510290.7  4157124;
            510325.0  4157133 ];
AOI.Zmin = 251.50;  AOI.Zmax = 276.93;

% Held-out sector: XY polygon and Z range. The upper limit (278.93 m) is
% above the top of the cloud (277.35 m): there is no upper limit.
polyHO = [ 510333.5  4157125;
           510291.1  4157115;
           510292.7  4157125;
           510325.0  4157133 ];
zHO = [268.50 278.93];

% Values of the paper, used for the checks
expect.nTLS = 65742;  expect.nHeld = 43776;  expect.nConstraint = 21966;
expect.nEval = 355;   expect.nCavity = 68;
expect.AUC_U = 0.614; expect.AUC_H = 0.614;  expect.AUC_C = 0.807;

%% 2) Read the input files
tlsFile = locate(root, {'Micio Conti_cropped_ERTarea.txt'});
fileU   = locate(root, {'Model_Micio Conti_U.txt'});
fileH   = locate(root, {'Model_Micio Conti_H_rho75000.txt'});
fileC   = locate(root, {'Model_Micio Conti_C_rho75000.txt'});

tls = readmatrix(tlsFile);  tls = tls(:,1:3);
if size(tls,1) ~= expect.nTLS
    error('quick_test:tls', 'The TLS cloud has %d points, %d expected: wrong file?', size(tls,1), expect.nTLS);
end
U = read_model(fileU);  H = read_model(fileH);  C = read_model(fileC);
if ~isequal(size(U), size(H)) || ~isequal(size(U), size(C)) || ...
        max(abs(U(:,1:3) - H(:,1:3)), [], 'all') > 1e-6 || max(abs(U(:,1:3) - C(:,1:3)), [], 'all') > 1e-6
    error('quick_test:mesh', 'U, H and C must have the same mesh.');
end

%% 3) Constraint region and held-out sector
inHO = inpolygon(tls(:,1), tls(:,2), polyHO(:,1), polyHO(:,2)) & tls(:,3) >= zHO(1) & tls(:,3) <= zHO(2);
heldPts = tls(inHO, :);
consPts = tls(~inHO, :);
lo = min(heldPts, [], 1);  hi = max(heldPts, [], 1);     % held-out rectangle
fprintf('TLS points: %d | constraint %d (%.1f%%) | held-out %d\n', ...
    size(tls,1), size(consPts,1), 100*size(consPts,1)/size(tls,1), size(heldPts,1));

%% 4) Evaluation set and ground truth
XYZ = U(:,1:3);
inTOT  = inpolygon(XYZ(:,1), XYZ(:,2), polyERT(:,1), polyERT(:,2)) & XYZ(:,3) >= AOI.Zmin & XYZ(:,3) <= AOI.Zmax;
inRect = all(XYZ >= lo & XYZ <= hi, 2);
sel = find(inTOT & inRect);
dC = minDistToSet(XYZ(sel,:), consPts);
dT = minDistToSet(XYZ(sel,:), tls);
sel = sel(dC > bufferRadius & dT <= maxDistTLS);
Xe = XYZ(sel,:);  n = numel(sel);

shp   = alphaShape(tls(:,1), tls(:,2), tls(:,3), alphaRadius);
truth = inShape(shp, Xe(:,1), Xe(:,2), Xe(:,3));
fprintf('Evaluation set: %d elements, %d inside the cavity\n', n, sum(truth));

%% 5) AUC on the evaluation set
rU = U(sel,4);  rH = H(sel,4);  rC = C(sel,4);
aucU = auc_rank(rU, truth);  aucH = auc_rank(rH, truth);  aucC = auc_rank(rC, truth);
fprintf('AUC  U = %.4f   H = %.4f   C = %.4f   (H-U = %+.4f, C-U = %+.4f)\n', ...
    aucU, aucH, aucC, aucH-aucU, aucC-aucU);

%% 6) Paired spatial block bootstrap
bc = floor((Xe - min(Xe,[],1)) ./ blockSize);
[~, ~, blockIndex] = unique(bc, 'rows');
nBlocks = max(blockIndex);
members = accumarray(blockIndex, (1:n)', [nBlocks,1], @(v){v}, {[]});
dH = nan(nBootstrap,1);  dCU = nan(nBootstrap,1);  k = 0;
while k < nBootstrap
    s   = randi(nBlocks, [nBlocks 1]);
    idx = cell2mat(members(s));
    t   = truth(idx);
    if all(t) || ~any(t), continue; end
    k = k + 1;
    a = auc_rank(rU(idx), t);
    dH(k)  = auc_rank(rH(idx), t) - a;
    dCU(k) = auc_rank(rC(idx), t) - a;
end
ciH = prctile(dH,  [2.5 97.5]);
ciC = prctile(dCU, [2.5 97.5]);
fprintf('Bootstrap (%d blocks, %d replicates)\n', nBlocks, nBootstrap);
fprintf('  Delta AUC H-U: median %+.3f  95%% interval [%.3f, %.3f]\n', median(dH),  ciH(1), ciH(2));
fprintf('  Delta AUC C-U: median %+.3f  95%% interval [%.3f, %.3f]\n', median(dCU), ciC(1), ciC(2));

%% 7) Comparison with the paper
checks = {
    'TLS points in the held-out sector',  size(heldPts,1) == expect.nHeld
    'TLS points in the constraint region', size(consPts,1) == expect.nConstraint
    'Evaluation elements',                 n == expect.nEval
    'Evaluation elements inside the cavity', sum(truth) == expect.nCavity
    'AUC of U',                            abs(aucU - expect.AUC_U) <= 0.002
    'AUC of H',                            abs(aucH - expect.AUC_H) <= 0.002
    'AUC of C',                            abs(aucC - expect.AUC_C) <= 0.005
    'Interval of H-U includes zero',       ciH(1) <= 0 && ciH(2) >= 0
    'Interval of C-U above zero',          ciC(1) > 0
    };
fprintf('\nChecks against the paper (Micio Conti, main evaluation set)\n');
allOK = true;
for i = 1:size(checks,1)
    if checks{i,2}, tag = 'PASS'; else, tag = 'FAIL'; allOK = false; end
    fprintf('  [%s] %s\n', tag, checks{i,1});
end
if allOK
    fprintf('\nAll checks passed in %.0f s.\n', toc(tStart));
else
    fprintf('\nSome checks failed. See the README, section "Quick test".\n');
end
% Expected bootstrap intervals in the paper (Table D.1, 3 m blocks, 100 blocks):
%   Delta AUC H-U  0.000 [-0.006, 0.006]     Delta AUC C-U  +0.192 [0.123, 0.269]
% The intervals can differ in the last digit because the random numbers differ.

%% ---------- local functions ----------
function f = locate(root, names)
% First existing file among 'names', searched anywhere under 'root'.
for i = 1:numel(names)
    d = dir(fullfile(root, '**', names{i}));
    d = d(~[d.isdir]);
    if ~isempty(d), f = fullfile(d(1).folder, d(1).name); return; end
end
error('quick_test:fileNotFound', 'File not found under %s: %s', root, strjoin(names, ' | '));
end

function M = read_model(f)
% Model file: X Y Z resistivity ... (with or without a header line).
M = readmatrix(f);
M = M(all(isfinite(M(:,1:4)), 2), 1:4);
end

function dmin = minDistToSet(P, Q)
% Minimum distance of each point of P from the set Q (chunked, no toolbox).
o = mean(Q, 1);  P = P - o;  Q = Q - o;  q2 = sum(Q.^2, 2)';
dmin = zeros(size(P,1), 1);  chunk = 200;
for i = 1:chunk:size(P,1)
    ii = i:min(i+chunk-1, size(P,1));  Pc = P(ii,:);
    D2 = sum(Pc.^2, 2) + q2 - 2*(Pc*Q');
    dmin(ii) = sqrt(max(min(D2, [], 2), 0));
end
end

function auc = auc_rank(score, truth)
% Exact ROC area = P(score of a positive > score of a negative), ties = 1/2.
nPos = sum(truth);  nNeg = sum(~truth);
if nPos == 0 || nNeg == 0, auc = NaN; return; end
[~, ~, g] = unique(score);  cnt = accumarray(g, 1);
avg = cumsum(cnt) - (cnt-1)/2;  ranks = avg(g);
auc = (sum(ranks(truth)) - nPos*(nPos+1)/2) / (nPos*nNeg);
end
