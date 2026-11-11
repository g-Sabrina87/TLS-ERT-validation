%% 6_2_kmeans_threshold_constrained.m
% Resistivity threshold that separates cavity and rock in an ERT model, with a
% K-means made robust in two ways:
%   1. TRIMMING: extreme outliers are excluded ONLY from the computation of the
%      centroids; every original element is then classified by the nearest
%      centroid on its TRUE value.
%   2. SPATIAL CONNECTIVITY FILTER: among the elements classified as cavity, only
%      the largest spatially connected block is kept.
% K-means is univariate (resistivity only). The threshold is then found by volume
% matching, against (a) the volume of the K-means cluster and (b) the exact volume
% of the TLS alpha shape.

clear; clc; close all;

%% ---------- 1) PARAMETERS ----------
modelFile = 'Model_Chiusazza_U.txt';        % ERT model (X Y Z resistivity)
siteName  = 'Chiusazza';

% Area investigated by the ERT
AOI.Xmin = 514191.9125;    AOI.Xmax = 514306.9915;
AOI.Ymin = 4097666.704;    AOI.Ymax = 4097785.755;
AOI.Zmin = 79.82597;  AOI.Zmax = 115.70548;

kRange = 1:10;

trimPercentile = 99;            % trimming percentile for the computation of the centroids
useConnectivityFilter = true;   % requires the Image Processing Toolbox (fallback included)

gridStep = 1.0;   % grid step for interpolation/isosurfaces (m)

% TLS geometry (for the non-circular comparison)
tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';   % full TLS cloud, same ERT area
alphaRadius = 2.5;

forceK = 3;                     % forced K ([] to use the K suggested by the silhouette)
includeHaloClusters = false;    % true = cavity = all clusters except the background

%% ---------- 2) READ AND CROP TO THE ERT AREA ----------
data = readmatrix(modelFile);
XYZ  = data(:,1:3);
rho  = data(:,4);

keep = XYZ(:,1) >= AOI.Xmin & XYZ(:,1) <= AOI.Xmax & ...
       XYZ(:,2) >= AOI.Ymin & XYZ(:,2) <= AOI.Ymax & ...
       XYZ(:,3) >= AOI.Zmin & XYZ(:,3) <= AOI.Zmax;
XYZ = XYZ(keep,:);
rho = rho(keep);

fprintf('Elements inside the ERT area: %d\n', numel(rho));

%% ---------- 3) ELBOW METHOD (inertia vs K) ----------
nK = numel(kRange);
inertia = nan(nK,1);

for i = 1:nK
    k = kRange(i);
    if k == 1
        inertia(i) = sum((rho - mean(rho)).^2);
    else
        [~, ~, sumd] = kmeans(rho, k, 'Replicates', 5);
        inertia(i) = sum(sumd);
    end
end

%% ---------- 4) SILHOUETTE SCORE (K >= 2) ----------
silScore = nan(nK,1);

maxPointsForSilhouette = 20000;
if numel(rho) > maxPointsForSilhouette
    rng(1);
    idxSample = randperm(numel(rho), maxPointsForSilhouette)';
    rhoForSil = rho(idxSample);
else
    rhoForSil = rho;
end

for i = 2:nK
    k = kRange(i);
    idxK = kmeans(rhoForSil, k, 'Replicates', 5);
    s = silhouette(rhoForSil, idxK);
    silScore(i) = mean(s);
end

%% ---------- 5) ELBOW + SILHOUETTE PLOTS ----------
figure;

subplot(2,1,1);
plot(kRange, inertia, 'o-', 'LineWidth', 1.5);
xlabel('Number of cluster (k)');
ylabel('Sum of squared distance');
title(sprintf('Elbow Method - %s', siteName));
grid on;

subplot(2,1,2);
plot(kRange(2:end), silScore(2:end), 'o-', 'LineWidth', 1.5);
xlabel('Number of cluster (k)');
ylabel('Silhouette Score');
title(sprintf('Silhouette Score - %s', siteName));
grid on;

[~, bestKidx] = max(silScore);
kOpt_silhouette = kRange(bestKidx);
fprintf('\nOptimal number of clusters (maximum Silhouette Score): K = %d\n', kOpt_silhouette);

if isempty(forceK)
    kOpt = kOpt_silhouette;
else
    kOpt = forceK;
    fprintf('K forced to %d (silhouette suggests %d).\n', kOpt, kOpt_silhouette);
end

%% ---------- 6) ROBUST K-MEANS WITH TRIMMING ----------
trimCutoff = prctile(rho, trimPercentile);
trimMask = rho <= trimCutoff;
rhoTrimmed = rho(trimMask);
fprintf('Trimming: %d elements above %.1f ohm*m (percentile %g) excluded ONLY from the computation of the centroids.\n', ...
    sum(~trimMask), trimCutoff, trimPercentile);

[~, centroids] = kmeans(rhoTrimmed, kOpt, 'Replicates', 10);

% Classify EVERY original element (including the outliers excluded above)
% by nearest centroid, on its TRUE value
sortedC = sort(centroids);
dists = abs(rho - sortedC');
[~, clusterIdx] = min(dists, [], 2);
centroids = sortedC;

[~, cavityClusterID] = max(centroids);   % highest-resistivity cluster = core

if includeHaloClusters && kOpt > 2
    isCavityCluster = clusterIdx ~= find(centroids == min(centroids), 1);
    fprintf('K=%d: CORE + HALO included as "cavity" (all clusters except the background).\n', kOpt);
else
    isCavityCluster = (clusterIdx == cavityClusterID);
    if kOpt > 2
        fprintf('K=%d: only the CORE (highest-resistivity cluster) is considered "cavity".\n', kOpt);
    end
end

fprintf('Cluster centroids (ohm*m):\n');
disp(centroids);
fprintf('Cavity cluster: #%d (centroid = %.1f ohm*m)\n', cavityClusterID, centroids(cavityClusterID));

%% ---------- 6b) CLUSTER SIZE ----------
fprintf('\n=== CLUSTER SIZE (after trimming) ===\n');
for c = 1:kOpt
    nInC = sum(clusterIdx == c);
    fprintf('Cluster #%d: %d elements (%.2f%% of the total), TRUE range [%.1f, %.1f] ohm*m\n', ...
        c, nInC, 100*nInC/numel(rho), min(rho(clusterIdx==c)), max(rho(clusterIdx==c)));
end

nCavity = sum(isCavityCluster);
if nCavity < 0.01 * numel(rho)
    warning('The cavity cluster contains less than 1%% of the elements even after trimming: try a lower trimPercentile.');
else
    fprintf('Cavity cluster size: %.2f%% of the elements.\n', 100*nCavity/numel(rho));
end

figure;
histogram(rho, 100);
hold on;
sortedCentroids = sort(centroids);
for i = 1:numel(sortedCentroids)-1
    xline(mean(sortedCentroids(i:i+1)), 'r--', 'LineWidth', 1.5);
end
xline(trimCutoff, 'g:', 'LineWidth', 1.5, 'Label', 'Trimming threshold');
xlabel('Resistivity (ohm\cdotm)');
ylabel('Number of elements');
set(gca, 'YScale', 'log');
title('Resistivity histogram with cluster boundary (red) and trimming threshold (green)');
grid on;

%% ---------- 7) INTERPOLATION ON A REGULAR GRID ----------
fprintf('\nInterpolation on a regular grid...\n');

xv = min(XYZ(:,1)):gridStep:max(XYZ(:,1));
yv = min(XYZ(:,2)):gridStep:max(XYZ(:,2));
zv = min(XYZ(:,3)):gridStep:max(XYZ(:,3));
[Xg, Yg, Zg] = meshgrid(xv, yv, zv);

Frho = scatteredInterpolant(XYZ(:,1), XYZ(:,2), XYZ(:,3), rho, 'linear', 'none');
Rhog = Frho(Xg, Yg, Zg);

Fclu = scatteredInterpolant(XYZ(:,1), XYZ(:,2), XYZ(:,3), double(isCavityCluster), 'linear', 'none');
Clug = Fclu(Xg, Yg, Zg);

fvCluster = isosurface(Xg, Yg, Zg, Clug, 0.5);
fprintf('Cluster boundary isosurface (before the connectivity filter): %d vertices\n', ...
    size(fvCluster.vertices,1));

if size(fvCluster.vertices,1) == 0
    error('Empty cluster isosurface: check gridStep or the quality of the clustering.');
end

%% ---------- 7b) SPATIAL CONNECTIVITY FILTER ----------
if useConnectivityFilter
    try
        clusterMask3D = Clug > 0.5;
        CC = bwconncomp(clusterMask3D, 26);

        compSizes = cellfun(@numel, CC.PixelIdxList);
        [maxSize, maxIdx] = max(compSizes);
        mainComponentIdx = CC.PixelIdxList{maxIdx};

        fprintf('\n=== SPATIAL CONNECTIVITY FILTER (bwconncomp) ===\n');
        fprintf('Number of connected components: %d\n', CC.NumObjects);
        fprintf('Component sizes (voxels, descending order):\n');
        disp(sort(compSizes, 'descend'));
    catch ME
        warning('bwconncomp not available or error (%s). Using the fallback without toolbox.', ME.message);

        clusterMask3D = Clug > 0.5;
        fprintf('\n=== SPATIAL CONNECTIVITY FILTER (fallback without toolbox) ===\n');
        [componentsList, compSizes] = conncomp3d_notoolbox(clusterMask3D);
        [maxSize, maxIdx] = max(compSizes);
        mainComponentIdx = componentsList{maxIdx};

        fprintf('Number of connected components: %d\n', numel(compSizes));
        fprintf('Component sizes (voxels, descending order):\n');
        disp(sort(compSizes, 'descend'));
    end

    fprintf('Main component: %d voxels out of %d in the cluster (%.1f%%)\n', ...
        maxSize, sum(clusterMask3D(:)), 100*maxSize/sum(clusterMask3D(:)));

    clusterMaskMain = false(size(clusterMask3D));
    clusterMaskMain(mainComponentIdx) = true;

    Clug = double(clusterMaskMain);
    fvCluster = isosurface(Xg, Yg, Zg, Clug, 0.5);

    fprintf('Cluster boundary isosurface (main component only): %d vertices\n', ...
        size(fvCluster.vertices,1));
end

%% ---------- 8) THRESHOLD BY VOLUME MATCHING ----------
lowBound  = min(centroids);
highBound = max(rho);
rhoRange = linspace(lowBound, highBound, 300);

cellVolume = gridStep^3;
clusterVolume = sum(Clug(:) > 0.5, 'omitnan') * cellVolume;
fprintf('\nVolume enclosed by the cluster (reference, after the connectivity filter): %.1f m^3\n', clusterVolume);

testVolume = nan(size(rhoRange));
matchDist  = nan(size(rhoRange));

fprintf('Searching the threshold with the best volume match...\n');
for i = 1:numel(rhoRange)
    testVolume(i) = sum(Rhog(:) >= rhoRange(i), 'omitnan') * cellVolume;

    fvTest = isosurface(Xg, Yg, Zg, Rhog, rhoRange(i));
    if size(fvTest.vertices,1) < 10
        continue;
    end
    d1 = pdist2(fvTest.vertices, fvCluster.vertices);
    matchDist(i) = mean([mean(min(d1,[],2)), mean(min(d1,[],1))]);
end

volDiff = abs(testVolume - clusterVolume);
[bestVolDiff, bestIdx] = min(volDiff);
thr_geomMatch = rhoRange(bestIdx);

fprintf('\n=== THRESHOLD BY VOLUME MATCHING (published values: 17000 or 4000 ohm*m) ===\n');
fprintf('Threshold: %.1f ohm*m\n', thr_geomMatch);
fprintf('Isosurface volume at this threshold: %.1f m^3 (vs cluster: %.1f m^3, difference %.1f%%)\n', ...
    testVolume(bestIdx), clusterVolume, 100*bestVolDiff/clusterVolume);
fprintf('Mean distance from the cluster boundary at this threshold: %.2f m\n', matchDist(bestIdx));

%% ---------- 8b) MATCHING CURVES ----------
figure;

subplot(1,2,1);
plot(rhoRange, testVolume, 'o-', 'LineWidth', 1.5);
hold on;
yline(clusterVolume, 'k--', 'Cluster volume (reference)');
plot(thr_geomMatch, testVolume(bestIdx), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
xlabel('Tested resistivity threshold (ohm\cdotm)');
ylabel('Enclosed volume (m^3)');
title('VOLUME matching');
legend('Isosurface volume', 'Cluster volume', 'Selected threshold', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(rhoRange, matchDist, 'o-', 'LineWidth', 1.5);
hold on;
plot(thr_geomMatch, matchDist(bestIdx), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
xlabel('Tested resistivity threshold (ohm\cdotm)');
ylabel('Mean distance from the cluster boundary (m)');
title('DISTANCE matching');
grid on;

%% ---------- 9) CLUSTER BOUNDARY VS ISOSURFACE ----------
fvBest = isosurface(Xg, Yg, Zg, Rhog, thr_geomMatch);

figure;
patch(fvCluster, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
patch(fvBest, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal; grid on; view(3);
camlight; lighting gouraud;
legend('K-means cluster boundary (after the connectivity filter)', ...
    sprintf('Isosurface at %.0f ohm\\cdotm', thr_geomMatch), 'Location', 'best');
title(sprintf('Cluster vs resistivity geometry - %s', siteName));
xlabel('X'); ylabel('Y'); zlabel('Z');

%% ============================================================
%  SECTION 10 - NON-CIRCULAR COMPARISON: threshold by direct matching with
%  the volume of the TLS alpha shape, instead of the K-means cluster.
%  ============================================================
fprintf('\n\n========== SECTION 10: COMPARISON WITH THE TLS VOLUME ==========\n');

%% ---------- 10a) TLS VOLUME (exact, from alphaShape) ----------
tlsData = readmatrix(tlsFullFile);
tlsXYZ  = tlsData(:,1:3);

shpTLS = alphaShape(tlsXYZ(:,1), tlsXYZ(:,2), tlsXYZ(:,3), alphaRadius);
volTLS_true = volume(shpTLS);

fprintf('TLS volume: %.1f m^3\n', volTLS_true);
fprintf('K-means cluster volume: %.1f m^3\n', clusterVolume);
fprintf('Cluster/TLS ratio: %.2fx\n', clusterVolume / volTLS_true);

%% ---------- 10b) THRESHOLD BY MATCHING WITH THE TLS VOLUME ----------
volDiff_TLS = abs(testVolume - volTLS_true);
[bestVolDiff_TLS, bestIdx_TLS] = min(volDiff_TLS);
thr_TLSMatch = rhoRange(bestIdx_TLS);

fvTLSMatch = isosurface(Xg, Yg, Zg, Rhog, thr_TLSMatch);
if size(fvTLSMatch.vertices,1) >= 10
    dTLS = pdist2(fvTLSMatch.vertices, tlsXYZ);
    matchDist_TLS = mean([mean(min(dTLS,[],2)), mean(min(dTLS,[],1))]);
else
    matchDist_TLS = NaN;
end

fprintf('\n=== THRESHOLD BY MATCHING WITH THE TLS VOLUME ===\n');
fprintf('Threshold: %.1f ohm*m\n', thr_TLSMatch);
fprintf('Isosurface volume at this threshold: %.1f m^3 (vs TLS: %.1f m^3, difference %.1f%%)\n', ...
    testVolume(bestIdx_TLS), volTLS_true, 100*bestVolDiff_TLS/volTLS_true);
fprintf('Mean distance from the TLS points: %.2f m\n', matchDist_TLS);

%% ---------- 10c) SUMMARY ----------
fprintf('\n=== SUMMARY (%s) ===\n', siteName);
fprintf('%-45s %12s\n', 'Method', 'Threshold (ohm*m)');
fprintf('%-45s %12.1f\n', 'Volume matching with the K-means cluster', thr_geomMatch);
fprintf('%-45s %12.1f\n', 'Volume matching with the TLS volume', thr_TLSMatch);
fprintf('%-45s %12s\n', 'Published value (Chiusazza)', '4000.0');

%% ---------- 10d) TLS, CLUSTER AND THE TWO THRESHOLDS ----------
figure;
plot(shpTLS, 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on;
patch(fvCluster, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch(fvBest, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
patch(fvTLSMatch, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
axis equal; grid on; view(3);
camlight; lighting gouraud;
legend('TLS', 'K-means cluster (core)', ...
    sprintf('Threshold from cluster (%.0f ohm\\cdotm)', thr_geomMatch), ...
    sprintf('Threshold from TLS (%.0f ohm\\cdotm)', thr_TLSMatch), ...
    'Location', 'best');
title(sprintf('TLS, cluster and thresholds - %s', siteName));
xlabel('X'); ylabel('Y'); zlabel('Z');

%% ---------- 10e) VOLUME vs THRESHOLD CURVE, WITH BOTH TARGETS ----------
figure;
plot(rhoRange, testVolume, 'o-', 'LineWidth', 1, 'MarkerSize', 3);
hold on;
yline(clusterVolume, 'b--', 'Cluster volume');
yline(volTLS_true, 'g--', 'TLS volume');
plot(thr_geomMatch, testVolume(bestIdx), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot(thr_TLSMatch, testVolume(bestIdx_TLS), 'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
xlabel('Tested resistivity threshold (ohm\cdotm)');
ylabel('Enclosed volume (m^3)');
title('Volume vs threshold');
legend('Isosurface volume', 'Cluster volume (target)', 'TLS volume (target)', ...
    'Selected threshold (cluster)', 'Selected threshold (TLS)', 'Location', 'best');
grid on;

%% ============================================================
%  LOCAL FUNCTION: 3D connected components without toolbox (fallback)
%  ============================================================
function [components, sizes] = conncomp3d_notoolbox(mask)
    [nx, ny, nz] = size(mask);
    visited = false(size(mask));
    idxList = find(mask);
    components = {};
    sizes = [];
    for ii = 1:numel(idxList)
        idx = idxList(ii);
        if visited(idx), continue; end
        queue = idx; visited(idx) = true; comp = idx; qi = 1;
        while qi <= numel(queue)
            cur = queue(qi); qi = qi + 1;
            [i,j,k] = ind2sub(size(mask), cur);
            for di = -1:1
                for dj = -1:1
                    for dk = -1:1
                        if di==0 && dj==0 && dk==0, continue; end
                        ni=i+di; nj=j+dj; nk=k+dk;
                        if ni>=1 && ni<=nx && nj>=1 && nj<=ny && nk>=1 && nk<=nz
                            nidx = sub2ind(size(mask), ni, nj, nk);
                            if mask(nidx) && ~visited(nidx)
                                visited(nidx) = true;
                                queue(end+1) = nidx; %#ok<AGROW>
                                comp(end+1) = nidx; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
        components{end+1} = comp; %#ok<AGROW>
        sizes(end+1) = numel(comp); %#ok<AGROW>
    end
end
