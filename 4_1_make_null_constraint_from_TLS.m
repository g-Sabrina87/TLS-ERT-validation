%% 4_1_make_null_constraint_from_TLS.m
% Generates NULL constraints: translated copies (no rotation by default) of the
% constraint point cloud, at random positions inside the ERT domain.
%
% Acceptance criteria:
%   - minimum distance of ANY point of the null from the full TLS cloud >= minDistanceFromReal
%     (Micio Conti: 3 m; Chiusazza: no constraint);
%   - fraction of points within nearRadius of the full TLS <= maxOverlapFraction;
%   - no point of the null inside the held-out rectangle enlarged by heldoutMargin;
%   - the null lies inside the mesh (no point farther than meshTol from an element centroid);
%   - vertical shift limited to +-dzMax, distance from the X/Y domain edges >= edgeMargin.
% Before generating, a feasibility scan estimates how many random positions are acceptable.
% startIndex allows generating replacements for discarded nulls without touching the existing ones.
%
% Output (outDir): TLS_NULL_NN.txt/.mat, NULL_summary_<site>.csv (shift, minimum distance
% from the TLS, 5th percentile, overlap), NULL_metadata_startNN.mat

clear; clc; close all;

%% ================= SETTINGS ==========================

site = 'Chiusazza';          % 'Micio Conti' or 'Chiusazza'
mode = 'reproduce';          % 'reproduce': rebuilds the nulls used in the paper from the stored translations
                             % 'generate' : draws new random translations (seed rngSeed)

switch site
    case 'Micio Conti'
        tlsFile       = 'Micio Conti_constraint_only.txt';   % cloud from which the nulls are generated
        realTlsFile   = 'Micio Conti_cropped_ERTarea.txt';   % FULL TLS cloud
        modelFile     = 'Model_Micio Conti_U.txt';           % gives ERT_DOMAIN and the mesh
        heldoutFile   = 'Micio Conti_holdout_region.mat';
        summaryFile   = 'NULL_summary_Micio Conti.csv';      % stored translations of the 30 nulls (second series)
        nNull         = 30;
        minDistanceFromReal = 6.0;   % m, no point of the null closer than this to the full TLS
        maxOverlapFraction  = 0.00;  % maximum fraction of points within nearRadius (0 = none)
        dzMax         = 4;           % m: maximum |vertical shift|
    case 'Chiusazza'
        tlsFile       = 'Chiusazza_constraint_only.txt';
        realTlsFile   = 'Chiusazza_cropped_ERTarea.txt';
        modelFile     = 'Model_Chiusazza_U.txt';
        heldoutFile   = '';                                  % no held-out check for the nulls of Chiusazza
        summaryFile   = 'NULL_summary_Chiusazza.csv';        % stored translations of the 20 nulls
        dzMax         = Inf;         % m: (the stored shifts reach 7.4 m vertically)
        nNull         = 20;
        minDistanceFromReal = 0;     % no constraint on the distance from the true cavity
        maxOverlapFraction  = 1.00;  % no limit on the overlap
end

outDir     = ['NULL_TLS_constraints_' strrep(site,' ','')];
startIndex = 1;           % first index of the file name
rngSeed    = 20261009;

nearRadius = 3.0;         % radius for the overlap fraction

% ERT domain [xmin xmax ymin ymax zmin zmax]: overwritten by the model file
ERT_DOMAIN = [0 0 0 0 0 0];

edgeMargin = 8;      % m: minimum distance of the null from the X/Y edges of the domain
meshTol   = 2.0;     % m
meshFrac  = 0.999;   % minimum fraction of points within meshTol

heldoutMargin = 1.5;    % m, enlargement of the held-out rectangle

feasibilityScan = true;
scanDistances   = [3 4.5 6];   % m
nScanTrials     = 1000;

maxTrialsPerNull = 50000;
nQuickCheck      = 2000;      % points of the quick check

useRandomRotation = false;    % translation only
maxRotationDeg    = 20;

writeXYZ = true;

%% ================================================================

if ~exist(outDir,'dir'), mkdir(outDir); end
rng(rngSeed,'twister');   % used only in 'generate' mode

%% 1. READ CLOUDS
fprintf('\nReading the constraint cloud: %s\n', tlsFile);
A = readmatrix(tlsFile);
if size(A,2) < 3, error('The file must have at least 3 columns (X Y Z).'); end
XYZ = A(:,1:3); XYZ = XYZ(all(isfinite(XYZ),2),:);

fprintf('Reading the full TLS cloud: %s\n', realTlsFile);
B = readmatrix(realTlsFile);
REAL = B(:,1:3); REAL = REAL(all(isfinite(REAL),2),:);

fprintf('Constraint points: %d | full TLS points: %d\n', size(XYZ,1), size(REAL,1));

tlsMin = min(XYZ,[],1); tlsMax = max(XYZ,[],1);
tlsCentroid = mean(XYZ,1); tlsSize = tlsMax - tlsMin;
fprintf('Size of the constraint cloud: %.2f x %.2f x %.2f m\n', tlsSize);

%% 1b. REPRODUCE MODE: rebuild the nulls from the stored translations
if strcmp(mode,'reproduce')
    S = readmatrix(summaryFile); S = S(all(isfinite(S(:,1:4)),2),1:4);   % columns: id, dx, dy, dz
    for k = 1:size(S,1)
        Pnull = XYZ + S(k,2:4);
        writematrix(Pnull, fullfile(outDir, sprintf('TLS_NULL_%02d.txt', S(k,1))), 'Delimiter','space');
        dAll = nnDist(Pnull, REAL);
        fprintf('NULL %02d | shift = [%7.2f %7.2f %7.2f] m | min distance from the TLS = %.2f m | overlap(%.1f m) = %.1f%%\n', ...
            S(k,1), S(k,2), S(k,3), S(k,4), min(dAll), nearRadius, 100*mean(dAll <= nearRadius));
    end
    fprintf('\nNulls rebuilt from %s into %s\n', summaryFile, outDir);
    return;
end

%% 2. DOMAIN (from the model) AND CONSISTENCY CHECK
MODEL = [];
if ~isempty(modelFile)
    Mx = readmatrix(modelFile); MODEL = Mx(:,1:3); MODEL = MODEL(all(isfinite(MODEL),2),:);
    ERT_DOMAIN = [min(MODEL(:,1)) max(MODEL(:,1)) min(MODEL(:,2)) max(MODEL(:,2)) ...
                  min(MODEL(:,3)) max(MODEL(:,3))];
    fprintf('ERT_DOMAIN from %s: X %.1f-%.1f, Y %.1f-%.1f, Z %.1f-%.1f\n', modelFile, ERT_DOMAIN);
end
if tlsMin(1) < ERT_DOMAIN(1)-5 || tlsMax(1) > ERT_DOMAIN(2)+5 || ...
   tlsMin(2) < ERT_DOMAIN(3)-5 || tlsMax(2) > ERT_DOMAIN(4)+5 || ...
   tlsMin(3) < ERT_DOMAIN(5)-5 || tlsMax(3) > ERT_DOMAIN(6)+5
    error(['The constraint cloud (X %.0f-%.0f, Y %.0f-%.0f, Z %.0f-%.0f) is NOT inside the ERT domain ', ...
           '(X %.0f-%.0f, Y %.0f-%.0f, Z %.0f-%.0f). Check modelFile / ERT_DOMAIN.'], ...
           tlsMin(1),tlsMax(1),tlsMin(2),tlsMax(2),tlsMin(3),tlsMax(3), ERT_DOMAIN);
end

domainMin  = [ERT_DOMAIN(1), ERT_DOMAIN(3), ERT_DOMAIN(5)];
domainMax  = [ERT_DOMAIN(2), ERT_DOMAIN(4), ERT_DOMAIN(6)];
domainSize = domainMax - domainMin;
if any(tlsSize > domainSize)
    error('The cloud is larger than the ERT domain: check ERT_DOMAIN.');
end

save(fullfile(outDir,'real_TLS_geometry.mat'), ...
    'XYZ','REAL','tlsMin','tlsMax','tlsCentroid','tlsSize','ERT_DOMAIN');

%% 2b. HELD-OUT AND FEASIBILITY SCAN
HO = [];
if ~isempty(heldoutFile) && isfile(heldoutFile)
    tmp = load(heldoutFile, 'HELDOUT'); HO = tmp.HELDOUT;
    HO.Xmin = HO.Xmin - heldoutMargin; HO.Xmax = HO.Xmax + heldoutMargin;
    HO.Ymin = HO.Ymin - heldoutMargin; HO.Ymax = HO.Ymax + heldoutMargin;
    HO.Zmin = HO.Zmin - heldoutMargin; HO.Zmax = HO.Zmax + heldoutMargin;
elseif ~isempty(heldoutFile)
    fprintf('WARNING: %s not found, held-out check skipped\n', heldoutFile);
end

if feasibilityScan
    fprintf('\nFeasibility scan (%d random translations):\n', nScanTrials);
    pMin0 = min(XYZ,[],1); pMax0 = max(XYZ,[],1);
    lowS = domainMin - pMin0; highS = domainMax - pMax0;
    lowS(1:2) = lowS(1:2) + edgeMargin; highS(1:2) = highS(1:2) - edgeMargin;
    lowS(3) = max(lowS(3), -dzMax); highS(3) = min(highS(3), dzMax);
    if any(lowS > highS), error('No room with the current dzMax/edgeMargin: reduce them.'); end
    dminScan = nan(nScanTrials,1); inHO = false(nScanTrials,1); inMesh = true(nScanTrials,1);
    for t = 1:nScanTrials
        sh = lowS + rand(1,3).*(highS - lowS);
        Pt = XYZ + sh;
        idq = randperm(size(Pt,1), min(nQuickCheck, size(Pt,1)));
        dminScan(t) = min(nnDist(Pt(idq,:), REAL));
        if ~isempty(MODEL), inMesh(t) = mean(nnDist(Pt(idq,:), MODEL) <= meshTol) >= meshFrac; end
        if ~isempty(HO)
            inHO(t) = any(Pt(:,1)>=HO.Xmin & Pt(:,1)<=HO.Xmax & Pt(:,2)>=HO.Ymin & ...
                          Pt(:,2)<=HO.Ymax & Pt(:,3)>=HO.Zmin & Pt(:,3)<=HO.Zmax);
        end
    end
    for dd = scanDistances
        ok = dminScan >= dd;
        fprintf('  minimum distance >= %.1f m: %.1f%% of the positions; also inside the mesh and outside the held-out sector: %.1f%%\n', ...
            dd, 100*mean(ok), 100*mean(ok & ~inHO & inMesh));
    end
    if mean(dminScan >= minDistanceFromReal & ~inHO & inMesh) < 0.001
        error('With minDistanceFromReal = %.1f m practically no position is valid.', minDistanceFromReal);
    end
end

%% 3. GENERATION
nullInfo = struct('id',{},'translation',{},'dminTrue',{},'p5True',{}, ...
                  'overlapFrac',{},'nPoints',{});
accepted = 0; trialGlobal = 0; nRejDist = 0; nRejOverlap = 0; nRejSize = 0; nRejMesh = 0;

fprintf('\nGenerating %d nulls (minimum distance >= %.2f m, overlap <= %.0f%% within %.2f m)...\n', ...
    nNull, minDistanceFromReal, 100*maxOverlapFraction, nearRadius);

while accepted < nNull
    trialGlobal = trialGlobal + 1;
    if trialGlobal > nNull*maxTrialsPerNull
        error(['Unable to generate %d valid nulls (rejected: %d by distance, ', ...
               '%d by overlap, %d by size).'], ...
               nNull, nRejDist, nRejOverlap, nRejSize);
    end

    P = XYZ; R = eye(3);
    if useRandomRotation
        ax = (2*rand-1)*maxRotationDeg*pi/180;
        ay = (2*rand-1)*maxRotationDeg*pi/180;
        az = (2*rand-1)*maxRotationDeg*pi/180;
        Rx = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
        Ry = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
        Rz = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
        R = Rz*Ry*Rx;
        P = (P - tlsCentroid)*R' + tlsCentroid;
    end

    pMin = min(P,[],1); pMax = max(P,[],1);
    if any((pMax - pMin) > domainSize), nRejSize = nRejSize + 1; continue; end

    % random translation such that ALL the points stay inside the domain
    lowShift  = domainMin - pMin;
    highShift = domainMax - pMax;
    lowShift(1:2)  = lowShift(1:2)  + edgeMargin;
    highShift(1:2) = highShift(1:2) - edgeMargin;
    lowShift(3)  = max(lowShift(3),  -dzMax);
    highShift(3) = min(highShift(3),  dzMax);
    if any(lowShift > highShift), nRejSize = nRejSize + 1; continue; end
    shift = lowShift + rand(1,3).*(highShift - lowShift);
    Pnull = P + shift;

    % --- quick check on a subset ---
    nQ  = min(nQuickCheck, size(Pnull,1));
    idx = randperm(size(Pnull,1), nQ);
    dq  = nnDist(Pnull(idx,:), REAL);
    if min(dq) < minDistanceFromReal, nRejDist = nRejDist + 1; continue; end

    % --- full check on ALL the points ---
    dAll = nnDist(Pnull, REAL);
    if min(dAll) < minDistanceFromReal, nRejDist = nRejDist + 1; continue; end
    overlapFrac = mean(dAll <= nearRadius);
    if overlapFrac > maxOverlapFraction, nRejOverlap = nRejOverlap + 1; continue; end
    if ~isempty(HO) && any(Pnull(:,1)>=HO.Xmin & Pnull(:,1)<=HO.Xmax & Pnull(:,2)>=HO.Ymin & ...
                           Pnull(:,2)<=HO.Ymax & Pnull(:,3)>=HO.Zmin & Pnull(:,3)<=HO.Zmax)
        nRejOverlap = nRejOverlap + 1; continue;
    end
    if ~isempty(MODEL) && mean(nnDist(Pnull, MODEL) <= meshTol) < meshFrac
        nRejMesh = nRejMesh + 1; continue;     % the null leaves the mesh
    end

    % --- accepted ---
    accepted = accepted + 1;
    fileIdx  = startIndex + accepted - 1;

    nullInfo(accepted).id          = fileIdx;
    nullInfo(accepted).translation = shift;
    nullInfo(accepted).dminTrue    = min(dAll);
    nullInfo(accepted).p5True      = prctile(dAll,5);
    nullInfo(accepted).overlapFrac = overlapFrac;
    nullInfo(accepted).nPoints     = size(Pnull,1);

    baseName = sprintf('TLS_NULL_%02d', fileIdx);
    if writeXYZ
        writematrix(Pnull, fullfile(outDir,[baseName '.txt']), 'Delimiter','space');
    end
    save(fullfile(outDir,[baseName '.mat']), 'Pnull','shift','R','ERT_DOMAIN');

    fprintf('NULL %02d | shift = [%7.2f %7.2f %7.2f] m | dmin = %.2f m | p5 = %.2f m | overlap(%.1f m) = %.1f%%\n', ...
        fileIdx, shift(1),shift(2),shift(3), min(dAll), prctile(dAll,5), nearRadius, 100*overlapFrac);
end

%% 4. SUMMARY
S = vertcat(nullInfo.translation);
T = table([nullInfo.id]', S(:,1), S(:,2), S(:,3), ...
          [nullInfo.dminTrue]', [nullInfo.p5True]', 100*[nullInfo.overlapFrac]', ...
          'VariableNames',{'NullID','dx_m','dy_m','dz_m', ...
                           'MinDistToTLS_m','P5DistToTLS_m','PctPointsWithinNearRadius'});
disp(T);
writetable(T, fullfile(outDir, sprintf('NULL_summary_%s.csv', strrep(site,' ',''))));

save(fullfile(outDir, sprintf('NULL_metadata_start%02d.mat', startIndex)), ...
    'nullInfo','nNull','rngSeed','minDistanceFromReal','nearRadius', ...
    'maxOverlapFraction','ERT_DOMAIN','useRandomRotation','maxRotationDeg', ...
    'tlsFile','realTlsFile','startIndex');

fprintf('\n=== DONE ===\nFolder: %s\nSeed: %d\nRejections: distance %d, overlap/held-out %d, outside mesh %d, size %d (out of %d trials)\n', ...
    outDir, rngSeed, nRejDist, nRejOverlap, nRejMesh, nRejSize, trialGlobal);

%% ================= LOCAL FUNCTION =============================
function d = nnDist(P, Q)
% Distance (m) of each point of P from the nearest point of Q.
    if exist('knnsearch','file') == 2
        [~, d] = knnsearch(Q, P, 'K', 1);   % the SECOND output is the distance
    else
        chunk = 1000;
        d = zeros(size(P,1),1);
        q2 = sum(Q.^2,2)';
        for i = 1:chunk:size(P,1)
            ii = i:min(i+chunk-1, size(P,1));
            Pc = P(ii,:);
            D2 = sum(Pc.^2,2) + q2 - 2*(Pc*Q');
            D2(D2 < 0) = 0;
            d(ii) = sqrt(min(D2,[],2));
        end
    end
end
