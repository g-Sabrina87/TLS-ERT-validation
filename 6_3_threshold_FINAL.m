%% 6_3_threshold_FINAL.m
% Resistivity threshold that delimits the cavity in an ERT model, with two methods:
%
%   METHOD 1 - thr_geomMatch: threshold whose isosurface has the same volume as the
%              K-means cluster (core, after outlier trimming + spatial connectivity filter)
%   METHOD 2 - thr_TLSMatch: threshold whose isosurface has the same volume as the
%              TLS geometry (independent, non-circular)

clear; clc; close all;

%% ============================================================
%  1) PARAMETERS
%  ============================================================

% ERT model to analyze (free or constrained): columns X Y Z Resistivity
modelFile = 'Model_Micio Conti_C_rho75000.txt';
siteName  = 'Micio Conti';

% Area investigated by the ERT (polygon + Z range)
polyXY = [ 510340.5  4157110;   % P1
           510311.1  4157106;   % P2
           510290.7  4157124;   % P3
           510325.0  4157133 ]; % P4

AOI.Zmin = 251.50;
AOI.Zmax = 276.93;

% TLS geometry (independent truth), cropped to the ERT polygon
tlsFullFile = 'Micio Conti_cropped_ERTarea.txt';
alphaRadius = 1.35;

% Value published in the original paper, for comparison
paperPublishedValue = 17000;   % ohm*m

kRange = 1:10;

% Preset: 'robusto' (K=3, trimming, connectivity filter) or 'originale'
% (K=2, no trimming, no filter: approach of the original paper)
analysisPreset = 'robusto';   % 'robusto' or 'originale'

switch analysisPreset
    case 'robusto'
        forceK = 3;
        includeHaloClusters = false;
        trimPercentile = 99;
        useConnectivityFilter = true;
    case 'originale'
        forceK = 2;
        includeHaloClusters = false;
        trimPercentile = 100;   % 100 = no trimming
        useConnectivityFilter = false;
    otherwise
        error('analysisPreset must be ''robusto'' or ''originale''.');
end
fprintf('Analysis preset: %s (K=%d, trimming=%g%%, connectivity=%d)\n', ...
    analysisPreset, forceK, trimPercentile, useConnectivityFilter);

gridStep = 1.0;   % m, grid for interpolation/isosurfaces

%% ============================================================
%  2) READ AND CROP TO THE ERT AREA
%  ============================================================
data = readmatrix(modelFile);
XYZ  = data(:,1:3);
rho  = data(:,4);

keep = inpolygon(XYZ(:,1), XYZ(:,2), polyXY(:,1), polyXY(:,2)) & ...
       XYZ(:,3) >= AOI.Zmin & XYZ(:,3) <= AOI.Zmax;
XYZ = XYZ(keep,:);
rho = rho(keep);

fprintf('=== %s ===\n', siteName);
fprintf('Total elements in the file (before cropping): %d\n', size(data,1));
fprintf('Elements inside the polygon + Z range (after cropping): %d\n', numel(rho));

if numel(rho) < 500
    warning('Very few elements inside the AOI (<500): check the limits.');
end

%% ============================================================
%  3) ELBOW METHOD + SILHOUETTE
%  ============================================================
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
    idxK = kmeans(rhoForSil, kRange(i), 'Replicates', 5);
    silScore(i) = mean(silhouette(rhoForSil, idxK));
end

figure('Name', 'Elbow + Silhouette');
subplot(2,1,1);
plot(kRange, inertia, 'o-', 'LineWidth', 1.5);
xlabel('Number of cluster (k)'); ylabel('Sum of squared distance');
title(sprintf('Elbow Method - %s', siteName)); grid on;

subplot(2,1,2);
plot(kRange(2:end), silScore(2:end), 'o-', 'LineWidth', 1.5);
xlabel('Number of cluster (k)'); ylabel('Silhouette Score');
title(sprintf('Silhouette Score - %s', siteName)); grid on;

[~, bestKidx] = max(silScore);
kOpt_silhouette = kRange(bestKidx);
kOpt = forceK;
fprintf('\nK suggested by the silhouette: %d | K used (forced): %d\n', kOpt_silhouette, kOpt);

%% ============================================================
%  4) ROBUST K-MEANS (OUTLIER TRIMMING) + IDENTIFICATION OF THE CORE
%  ============================================================
trimCutoff = prctile(rho, trimPercentile);
trimMask = rho <= trimCutoff;
rhoTrimmed = rho(trimMask);
fprintf('Trimming: %d elements above %.1f ohm*m (percentile %g) excluded from the computation of the centroids.\n', ...
    sum(~trimMask), trimCutoff, trimPercentile);

[~, centroids] = kmeans(rhoTrimmed, kOpt, 'Replicates', 10);
sortedC = sort(centroids);
dists = abs(rho - sortedC');
[~, clusterIdx] = min(dists, [], 2);
centroids = sortedC;

[~, cavityClusterID] = max(centroids);
if includeHaloClusters && kOpt > 2
    isCavityCluster = clusterIdx ~= find(centroids == min(centroids), 1);
else
    isCavityCluster = (clusterIdx == cavityClusterID);
end

fprintf('\nCluster centroids (ohm*m):\n'); disp(centroids);
fprintf('Cavity cluster: #%d (centroid = %.1f ohm*m)\n', cavityClusterID, centroids(cavityClusterID));

fprintf('\n--- Cluster size ---\n');
for c = 1:kOpt
    nInC = sum(clusterIdx == c);
    fprintf('Cluster #%d: %d elements (%.2f%%), range [%.1f, %.1f] ohm*m\n', ...
        c, nInC, 100*nInC/numel(rho), min(rho(clusterIdx==c)), max(rho(clusterIdx==c)));
end

figure('Name', 'Resistivity histogram');
histogram(rho, 100);
hold on;
for i = 1:numel(centroids)-1
    xline(mean(centroids(i:i+1)), 'r--', 'LineWidth', 1.5);
end
xline(trimCutoff, 'g:', 'LineWidth', 1.5, 'Label', 'Trimming threshold');
set(gca, 'YScale', 'log');
xlabel('Resistivity (ohm\cdotm)'); ylabel('Number of elements');
title(sprintf('Resistivity histogram - %s', siteName)); grid on;

%% ============================================================
%  5) INTERPOLATION ON A REGULAR GRID
%  ============================================================
xv = min(XYZ(:,1)):gridStep:max(XYZ(:,1));
yv = min(XYZ(:,2)):gridStep:max(XYZ(:,2));
zv = min(XYZ(:,3)):gridStep:max(XYZ(:,3));
[Xg, Yg, Zg] = meshgrid(xv, yv, zv);

Frho = scatteredInterpolant(XYZ(:,1), XYZ(:,2), XYZ(:,3), rho, 'linear', 'none');
Rhog = Frho(Xg, Yg, Zg);

Fclu = scatteredInterpolant(XYZ(:,1), XYZ(:,2), XYZ(:,3), double(isCavityCluster), 'linear', 'none');
Clug = Fclu(Xg, Yg, Zg);

fvCluster = isosurface(Xg, Yg, Zg, Clug, 0.5);
if isempty(fvCluster.vertices)
    error('Empty cluster isosurface: check gridStep or the quality of the clustering.');
end

%% ============================================================
%  6) SPATIAL CONNECTIVITY FILTER (keeps only the compact block)
%  ============================================================
if useConnectivityFilter
    clusterMask3D = Clug > 0.5;
    try
        CC = bwconncomp(clusterMask3D, 26);
        compSizes = cellfun(@numel, CC.PixelIdxList);
        [maxSize, maxIdx] = max(compSizes);
        mainComponentIdx = CC.PixelIdxList{maxIdx};
        methodUsed = 'bwconncomp';
    catch ME
        warning('bwconncomp not available (%s). Using the fallback without toolbox.', ME.message);
        [componentsList, compSizes] = conncomp3d_notoolbox(clusterMask3D);
        [maxSize, maxIdx] = max(compSizes);
        mainComponentIdx = componentsList{maxIdx};
        methodUsed = 'fallback';
    end

    fprintf('\n--- Connectivity filter (%s) ---\n', methodUsed);
    fprintf('Components found: %d | Main component: %d/%d voxels (%.1f%%)\n', ...
        numel(compSizes), maxSize, sum(clusterMask3D(:)), 100*maxSize/sum(clusterMask3D(:)));

    clusterMaskMain = false(size(clusterMask3D));
    clusterMaskMain(mainComponentIdx) = true;
    Clug = double(clusterMaskMain);
    fvCluster = isosurface(Xg, Yg, Zg, Clug, 0.5);
end

cellVolume = gridStep^3;
clusterVolume = sum(Clug(:) > 0.5, 'omitnan') * cellVolume;
fprintf('\nCluster volume (METHOD 1 reference): %.1f m^3\n', clusterVolume);

%% ============================================================
%  7) TLS VOLUME (independent truth, for METHOD 2)
%  ============================================================
tlsData = readmatrix(tlsFullFile);
tlsXYZ  = tlsData(:,1:3);
shpTLS  = alphaShape(tlsXYZ(:,1), tlsXYZ(:,2), tlsXYZ(:,3), alphaRadius);
volTLS_true = volume(shpTLS);
fprintf('TLS volume (METHOD 2 reference): %.1f m^3\n', volTLS_true);
fprintf('Cluster/TLS ratio: %.2fx\n', clusterVolume/volTLS_true);

%% ============================================================
%  8) SEARCH OF THE TWO THRESHOLDS
%  ============================================================
lowBound  = prctile(rho, 1);
highBound = prctile(rho, 99.9);
rhoRange  = linspace(lowBound, highBound, 300);

testVolume = nan(size(rhoRange));
for i = 1:numel(rhoRange)
    testVolume(i) = sum(Rhog(:) >= rhoRange(i), 'omitnan') * cellVolume;
end

% METHOD 1: matching with the cluster
[diff1, idx1] = min(abs(testVolume - clusterVolume));
thr_geomMatch = rhoRange(idx1);

% METHOD 2: matching with the TLS volume
[diff2, idx2] = min(abs(testVolume - volTLS_true));
thr_TLSMatch = rhoRange(idx2);

fvBest    = isosurface(Xg, Yg, Zg, Rhog, thr_geomMatch);
fvTLSBest = isosurface(Xg, Yg, Zg, Rhog, thr_TLSMatch);

fprintf('\n========== FINAL RESULTS - %s ==========\n', siteName);
fprintf('%-48s %14s %14s\n', 'Method', 'Threshold(ohm*m)', 'Vol. diff.(%)');
fprintf('%-48s %14.1f %14.1f\n', 'METHOD 1 - matching with the K-means cluster', ...
    thr_geomMatch, 100*diff1/clusterVolume);
fprintf('%-48s %14.1f %14.1f\n', 'METHOD 2 - matching with the TLS volume', ...
    thr_TLSMatch, 100*diff2/volTLS_true);
fprintf('%-48s %14.1f %14s\n', 'Value published in the original paper', paperPublishedValue, '-');

%% ============================================================
%  9) SUMMARY PLOTS
%  ============================================================
figure('Name', 'Final comparison');
plot(shpTLS, 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
patch(fvCluster, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch(fvBest, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
patch(fvTLSBest, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
axis equal; grid on; view(3); camlight; lighting gouraud;
legend('TLS', 'K-means cluster (core)', ...
    sprintf('METHOD 1 threshold (%.0f ohm\\cdotm)', thr_geomMatch), ...
    sprintf('METHOD 2 threshold (%.0f ohm\\cdotm)', thr_TLSMatch), ...
    'Location', 'best');
title(sprintf('%s - final comparison', siteName));
xlabel('X'); ylabel('Y'); zlabel('Z');

figure('Name', 'Volume vs threshold');
plot(rhoRange, testVolume, 'o-', 'MarkerSize', 3); hold on;
yline(clusterVolume, 'b--', 'Cluster volume');
yline(volTLS_true, 'g--', 'TLS volume');
plot(thr_geomMatch, testVolume(idx1), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot(thr_TLSMatch, testVolume(idx2), 'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
xline(paperPublishedValue, 'k:', 'LineWidth', 1.5, 'Label', 'Published value');
xlabel('Threshold (ohm\cdotm)'); ylabel('Enclosed volume (m^3)');
title('Volume vs threshold curve'); grid on;
legend('Isosurface volume', 'Cluster target', 'TLS target', ...
    'METHOD 1 threshold', 'METHOD 2 threshold', 'Location', 'best');

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
