%% 6_4_threshold_sensitivity.m
% Sensitivity analysis: repeats the threshold pipeline (Method 1 - K-means cluster;
% Method 2 - TLS volume) for several constrained models, each with a different
% resistivity imposed on the cavity anomaly (rho_input), and produces a table and a
% threshold-vs-rho_input plot.

clear; clc; close all;

%% ============================================================
%  1) PARAMETERS
%  ============================================================

% Runs to compare: one file per imposed rho value
runs = struct('file', {}, 'rho_input', {});
runs(1).file = 'Model_Micio Conti_C_rho17000.txt';   runs(1).rho_input = 17000;
runs(2).file = 'Model_Micio Conti_C_rho50000.txt';   runs(2).rho_input = 50000;
runs(3).file = 'Model_Micio Conti_C_rho75000.txt';   runs(3).rho_input = 75000;
runs(4).file = 'Model_Micio Conti_C_rho100000.txt';  runs(4).rho_input = 100000;

siteName = 'Micio Conti';

% Area investigated by the ERT (polygon + Z range)
polyXY = [ 510340.5  4157110;   % P1
           510311.1  4157106;   % P2
           510290.7  4157124;   % P3
           510325.0  4157133 ]; % P4
AOI.Zmin = 251.50;
AOI.Zmax = 276.93;

% TLS geometry (independent truth, computed once), cropped to the ERT polygon
tlsFullFile = 'Micio Conti_cropped_ERTarea.txt';
alphaRadius = 1.35;

paperPublishedValue = 17000;   % ohm*m, reference in the plots

% K-means preset
% 'robusto'   -> K=3, trimming 99%, connectivity filter
% 'originale' -> K=2, no trimming, no filter (approach of the original paper)
analysisPreset = 'originale';

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
fprintf('Analysis preset: %s (K=%d, trimming=%g%%, connectivity=%d)\n\n', ...
    analysisPreset, forceK, trimPercentile, useConnectivityFilter);

gridStep = 1.0;

%% ============================================================
%  2) TLS VOLUME (computed once, same for all the runs)
%  ============================================================
tlsData = readmatrix(tlsFullFile);
tlsXYZ  = tlsData(:,1:3);
shpTLS  = alphaShape(tlsXYZ(:,1), tlsXYZ(:,2), tlsXYZ(:,3), alphaRadius);
volTLS_true = volume(shpTLS);
fprintf('TLS volume (common to all the runs): %.1f m^3\n\n', volTLS_true);

%% ============================================================
%  3) LOOP OVER THE RUNS
%  ============================================================
nRuns = numel(runs);
results = table('Size', [nRuns 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'rho_input','clusterVolume','thr_geomMatch','thr_TLSMatch','ratio_cluster_TLS'});

for r = 1:nRuns
    fprintf('\n========== RUN %d/%d: rho_input = %g ohm-m ==========\n', r, nRuns, runs(r).rho_input);

    % --- Read and crop to the polygon ---
    data = readmatrix(runs(r).file);
    XYZ  = data(:,1:3);
    rho  = data(:,4);

    keep = inpolygon(XYZ(:,1), XYZ(:,2), polyXY(:,1), polyXY(:,2)) & ...
           XYZ(:,3) >= AOI.Zmin & XYZ(:,3) <= AOI.Zmax;
    XYZ = XYZ(keep,:);
    rho = rho(keep);
    fprintf('Elements inside the ERT area: %d\n', numel(rho));

    % --- Robust K-means (trimming) ---
    trimCutoff = prctile(rho, trimPercentile);
    rhoTrimmed = rho(rho <= trimCutoff);
    [~, centroids] = kmeans(rhoTrimmed, forceK, 'Replicates', 10);
    sortedC = sort(centroids);
    [~, clusterIdx] = min(abs(rho - sortedC'), [], 2);
    centroids = sortedC;

    [~, cavityClusterID] = max(centroids);
    if includeHaloClusters && forceK > 2
        isCavityCluster = clusterIdx ~= find(centroids == min(centroids), 1);
    else
        isCavityCluster = (clusterIdx == cavityClusterID);
    end
    fprintf('Cavity cluster: #%d, centroid = %.1f ohm*m, %d elements (%.2f%%)\n', ...
        cavityClusterID, centroids(cavityClusterID), sum(isCavityCluster), 100*sum(isCavityCluster)/numel(rho));

    % --- Interpolation on a grid ---
    xv = min(XYZ(:,1)):gridStep:max(XYZ(:,1));
    yv = min(XYZ(:,2)):gridStep:max(XYZ(:,2));
    zv = min(XYZ(:,3)):gridStep:max(XYZ(:,3));
    [Xg, Yg, Zg] = meshgrid(xv, yv, zv);

    Frho = scatteredInterpolant(XYZ(:,1), XYZ(:,2), XYZ(:,3), rho, 'linear', 'none');
    Rhog = Frho(Xg, Yg, Zg);

    Fclu = scatteredInterpolant(XYZ(:,1), XYZ(:,2), XYZ(:,3), double(isCavityCluster), 'linear', 'none');
    Clug = Fclu(Xg, Yg, Zg);

    % --- Connectivity filter ---
    if useConnectivityFilter
        clusterMask3D = Clug > 0.5;
        try
            CC = bwconncomp(clusterMask3D, 26);
            compSizes = cellfun(@numel, CC.PixelIdxList);
            [~, maxIdx] = max(compSizes);
            mainComponentIdx = CC.PixelIdxList{maxIdx};
        catch
            [componentsList, compSizes] = conncomp3d_notoolbox(clusterMask3D);
            [~, maxIdx] = max(compSizes);
            mainComponentIdx = componentsList{maxIdx};
        end
        clusterMaskMain = false(size(clusterMask3D));
        clusterMaskMain(mainComponentIdx) = true;
        Clug = double(clusterMaskMain);
    end

    cellVolume = gridStep^3;
    clusterVolume = sum(Clug(:) > 0.5, 'omitnan') * cellVolume;

    % --- Search of the two thresholds ---
    lowBound  = prctile(rho, 1);
    highBound = prctile(rho, 99.9);
    rhoRange  = linspace(lowBound, highBound, 300);

    testVolume = nan(size(rhoRange));
    for i = 1:numel(rhoRange)
        testVolume(i) = sum(Rhog(:) >= rhoRange(i), 'omitnan') * cellVolume;
    end

    [~, idx1] = min(abs(testVolume - clusterVolume));
    thr_geomMatch = rhoRange(idx1);

    [~, idx2] = min(abs(testVolume - volTLS_true));
    thr_TLSMatch = rhoRange(idx2);

    fprintf('Cluster volume: %.1f m^3 | thr_geomMatch = %.1f ohm-m | thr_TLSMatch = %.1f ohm-m\n', ...
        clusterVolume, thr_geomMatch, thr_TLSMatch);

    results.rho_input(r)        = runs(r).rho_input;
    results.clusterVolume(r)    = clusterVolume;
    results.thr_geomMatch(r)    = thr_geomMatch;
    results.thr_TLSMatch(r)     = thr_TLSMatch;
    results.ratio_cluster_TLS(r)= clusterVolume / volTLS_true;
end

%% ============================================================
%  4) SUMMARY AND SENSITIVITY PLOT
%  ============================================================
fprintf('\n\n========== SENSITIVITY SUMMARY - %s ==========\n', siteName);
disp(results);

figure('Name', 'Threshold sensitivity vs imposed rho');
plot(results.rho_input, results.thr_geomMatch, 'bs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
hold on;
plot(results.rho_input, results.thr_TLSMatch, 'g^-', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'MarkerSize', 8);
yline(paperPublishedValue, 'k--', 'LineWidth', 1.5, 'Label', sprintf('Published (%.0f)', paperPublishedValue));
set(gca, 'XScale', 'log');
xlabel('Resistivity imposed on the constraint, rho\_input (ohm\cdotm)');
ylabel('Threshold (ohm\cdotm)');
title(sprintf('Sensitivity of the threshold to the imposed constraint value - %s', siteName));
legend('Method 1 (K-means cluster)', 'Method 2 (TLS volume)', 'Published value', 'Location', 'best');
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
