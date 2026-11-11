%% 2_2_J_alpha_weight_formalization.m
% J(alpha) with explicit normalized weights:
%
%     J(alpha) = wD*D*(alpha) + wT*T*(alpha) + wC*C*(alpha),   wD+wT+wC = 1
%     each term min-max normalized in [0,1]: D* = (D-Dmin)/(Dmax-Dmin), same for T*, C*
%
% 1. Raw components for a range of alphaRadius, at both sites:
%    D = fidelity (Chamfer), T = topological penalty, C = voxel-coverage fragmentation.
% 2. Min-max normalization and J with equal weights (1/3, 1/3, 1/3).
% 3. Sensitivity of alpha* = argmin J to the weights (each weight in [0.2, 0.6], sum = 1).
% 4. Hierarchical criterion: only alphas with T = 0, then minimize wD*D* + wC*C*
%    (renormalized on the valid subset), with sensitivity to wD.

clear; clc; close all;

%% ============================================================
%  1) SITE CONFIGURATION
%  ============================================================

sites = struct('name', {}, 'tlsFile', {}, 'polyXY', {}, 'AOIZmin', {}, 'AOIZmax', {}, ...
    'alphaTest', {}, 'voxelSize', {});

% --- Micio Conti ---
sites(1).name = 'Micio Conti';
sites(1).tlsFile = 'Micio Conti_cropped_ERTarea.txt';
sites(1).polyXY = [ 510340.5  4157110;
                     510311.1  4157106;
                     510290.7  4157124;
                     510325.0  4157133 ];
sites(1).AOIZmin = 251.50;
sites(1).AOIZmax = 276.93;
sites(1).alphaTest = 0.75:0.1:2.0;
sites(1).voxelSize = 0.5;

% --- Chiusazza ---
sites(2).name = 'Chiusazza';
sites(2).tlsFile = 'Chiusazza_cropped_ERTarea.txt';
sites(2).polyXY = [ 514191.91  4097666.70;
                     514306.99  4097666.70;
                     514306.99  4097785.75;
                     514191.91  4097785.75 ];
sites(2).AOIZmin = 79.82;
sites(2).AOIZmax = 115.71;
sites(2).alphaTest = 0.5:0.5:6.0;
sites(2).voxelSize = 1.0;

maxPointsForFidelity = 5000;

%% ============================================================
%  2) RAW COMPONENTS FOR EACH SITE
%  ============================================================

allResults = struct();

for s = 1:numel(sites)
    site = sites(s);
    fprintf('\n\n========== %s: raw components of J(alpha) ==========\n', site.name);

    data = readmatrix(site.tlsFile);
    xyz  = data(:,1:3);
    fprintf('Points read: %d\n', size(xyz,1));

    capBounds = [min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) min(xyz(:,3)) max(xyz(:,3))];
    polyBounds = [min(site.polyXY(:,1)) max(site.polyXY(:,1)) min(site.polyXY(:,2)) max(site.polyXY(:,2))];
    box0 = [max(polyBounds(1), capBounds(1)), min(polyBounds(2), capBounds(2)), ...
            max(polyBounds(3), capBounds(3)), min(polyBounds(4), capBounds(4)), ...
            max(site.AOIZmin, capBounds(5)), min(site.AOIZmax, capBounds(6))];

    xv = box0(1):site.voxelSize:box0(2);
    yv = box0(3):site.voxelSize:box0(4);
    zv = box0(5):site.voxelSize:box0(6);
    [XX, YY, ZZ] = ndgrid(xv, yv, zv);
    insidePolyMask = inpolygon(XX(:), YY(:), site.polyXY(:,1), site.polyXY(:,2));
    insidePolyMask = reshape(insidePolyMask, size(XX));

    if size(xyz,1) > maxPointsForFidelity
        rng(1);
        idxSample = randperm(size(xyz,1), maxPointsForFidelity)';
        xyzForFidelity = xyz(idxSample,:);
    else
        xyzForFidelity = xyz;
    end

    nTest = numel(site.alphaTest);
    D_raw = nan(nTest,1);   % fidelity (Chamfer)
    T_raw = nan(nTest,1);   % raw topological penalty (regions - 1)
    C_raw = nan(nTest,1);   % voxel-coverage fragmentation (fraction)
    volTest = nan(nTest,1);

    for i = 1:nTest
        shpTest = alphaShape(xyz(:,1), xyz(:,2), xyz(:,3), site.alphaTest(i));
        volTest(i) = volume(shpTest);
        try
            nR = numRegions(shpTest);
        catch
            nR = NaN;
        end

        if volTest(i) < 1e-6 || isnan(nR)
            fprintf('  alpha=%.2f -> degenerate shape, discarded\n', site.alphaTest(i));
            continue;
        end
        T_raw(i) = max(0, nR - 1);

        try
            [tri, triPts] = boundaryFacets(shpTest);
        catch
            tri = [];
        end
        if isempty(tri)
            continue;
        end
        nSamplesPerTri = 3;
        v1 = triPts(tri(:,1),:); v2 = triPts(tri(:,2),:); v3 = triPts(tri(:,3),:);
        surfSamples = zeros(size(tri,1)*nSamplesPerTri, 3);
        k = 0;
        for t = 1:size(tri,1)
            for ss = 1:nSamplesPerTri
                r1 = rand(); r2 = rand();
                if r1 + r2 > 1, r1 = 1-r1; r2 = 1-r2; end
                k = k + 1;
                surfSamples(k,:) = v1(t,:) + r1*(v2(t,:)-v1(t,:)) + r2*(v3(t,:)-v1(t,:));
            end
        end
        surfSamples = surfSamples(1:k,:);
        if size(surfSamples,1) > maxPointsForFidelity
            idxS = randperm(size(surfSamples,1), maxPointsForFidelity)';
            surfSamples = surfSamples(idxS,:);
        end
        d1 = pairwiseDistances_notoolbox(xyzForFidelity, surfSamples);
        D_raw(i) = mean([mean(min(d1,[],2)), mean(min(d1,[],1))]);

        insideCavity = inShape(shpTest, XX(:), YY(:), ZZ(:));
        insideCavity = reshape(insideCavity, size(XX));
        occMask = insideCavity & insidePolyMask;
        nOccupied = nnz(occMask);
        if nOccupied == 0
            C_raw(i) = 1;
        else
            try
                CC = bwconncomp(occMask, 26);
                compSizes = cellfun(@numel, CC.PixelIdxList);
            catch
                [~, compSizes] = conncomp3d_notoolbox(occMask);
            end
            mainSize = max(compSizes);
            C_raw(i) = 1 - mainSize/nOccupied;
        end

        fprintf('  alpha=%.2f -> D=%.3f m, T=%d (extra regions), C=%.1f%%\n', ...
            site.alphaTest(i), D_raw(i), T_raw(i), 100*C_raw(i));
    end

    allResults.(matlab.lang.makeValidName(site.name)).alphaTest = site.alphaTest;
    allResults.(matlab.lang.makeValidName(site.name)).D_raw = D_raw;
    allResults.(matlab.lang.makeValidName(site.name)).T_raw = T_raw;
    allResults.(matlab.lang.makeValidName(site.name)).C_raw = C_raw;
end

%% ============================================================
%  3) NORMALIZATION, J WITH EQUAL WEIGHTS, WEIGHT SENSITIVITY
%  ============================================================

weightStep = 0.05;
weightMin = 0.2;
weightMax = 0.6;

figure('Name', 'J(alpha) - equal weights');
figure('Name', 'Weight sensitivity - alpha*');
figure('Name', 'D/C sensitivity - hierarchical criterion');

for s = 1:numel(sites)
    siteName = sites(s).name;
    fieldName = matlab.lang.makeValidName(siteName);
    R = allResults.(fieldName);
    alphaTest = R.alphaTest;

    valid = ~isnan(R.D_raw) & ~isnan(R.T_raw) & ~isnan(R.C_raw);
    alphaValid = alphaTest(valid);
    D = R.D_raw(valid); T = R.T_raw(valid); C = R.C_raw(valid);

    % --- Min-max normalization ---
    Dstar = (D - min(D)) / max(range(D), eps);
    Tstar = (T - min(T)) / max(range(T), eps);
    Cstar = (C - min(C)) / max(range(C), eps);

    % --- J with equal weights ---
    wD = 1/3; wT = 1/3; wC = 1/3;
    J_equal = wD*Dstar + wT*Tstar + wC*Cstar;
    [Jmin, idxMin] = min(J_equal);
    alphaOptimal_equal = alphaValid(idxMin);

    fprintf('\n=== %s: J(alpha) with equal weights (1/3,1/3,1/3) ===\n', siteName);
    fprintf('alpha* = %.2f (J=%.3f)\n', alphaOptimal_equal, Jmin);

    figure(findobj('Name', 'J(alpha) - equal weights'));
    subplot(1, numel(sites), s);
    plot(alphaValid, J_equal, 'o-', 'LineWidth', 1.5); hold on;
    plot(alphaOptimal_equal, Jmin, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
    xlabel('alphaRadius (m)'); ylabel('J(alpha), equal weights');
    title(sprintf('%s (alpha*=%.2f)', siteName, alphaOptimal_equal));
    grid on;

    % --- Weight sensitivity ---
    wGrid = weightMin:weightStep:weightMax;
    alphaStarList = [];
    weightCombos = [];

    for wDi = wGrid
        for wTi = wGrid
            wCi = 1 - wDi - wTi;
            if wCi < weightMin || wCi > weightMax
                continue;   % combination outside the allowed range
            end
            Jw = wDi*Dstar + wTi*Tstar + wCi*Cstar;
            [~, idxW] = min(Jw);
            alphaStarList(end+1) = alphaValid(idxW); %#ok<AGROW>
            weightCombos(end+1,:) = [wDi, wTi, wCi]; %#ok<AGROW>
        end
    end

    nCombos = numel(alphaStarList);
    modeAlphaStar = mode(alphaStarList);
    fracStable = mean(alphaStarList == modeAlphaStar);
    rangeAlphaStar = [min(alphaStarList), max(alphaStarList)];

    fprintf('Weight sensitivity (%d valid combinations, each weight in [%.1f, %.1f]):\n', ...
        nCombos, weightMin, weightMax);
    fprintf('  most frequent alpha*: %.2f (in %.1f%% of the combinations)\n', modeAlphaStar, 100*fracStable);
    fprintf('  observed range of alpha*: [%.2f, %.2f]\n', rangeAlphaStar(1), rangeAlphaStar(2));

    figure(findobj('Name', 'Weight sensitivity - alpha*'));
    subplot(1, numel(sites), s);
    histogram(alphaStarList, 'BinMethod', 'integers');
    hold on;
    xline(alphaOptimal_equal, 'r--', 'LineWidth', 2, 'Label', 'Equal weights');
    xlabel('alpha*'); ylabel('Number of weight combinations');
    title(sprintf('%s: stability of alpha* (%.0f%% at the mode)', siteName, 100*fracStable));
    grid on;

    allResults.(fieldName).Dstar = Dstar;
    allResults.(fieldName).Tstar = Tstar;
    allResults.(fieldName).Cstar = Cstar;
    allResults.(fieldName).alphaOptimal_equal = alphaOptimal_equal;
    allResults.(fieldName).alphaStarList = alphaStarList;
    allResults.(fieldName).fracStable = fracStable;

    %% ---------- HIERARCHICAL CRITERION ----------
    % Step 1: keep only the topologically valid alphas (T = 0, not normalized).
    % Step 2: among these, minimize the weighted combination of D and C
    % (renormalized on the valid subset only).
    fprintf('\n--- %s: HIERARCHICAL CRITERION (T = 0 as prerequisite) ---\n', siteName);

    validTopo = (T == 0);
    if ~any(validTopo)
        warning('%s: no tested alpha has T = 0 - extend the range of alphaTest.', siteName);
    else
        alphaValidTopo = alphaValid(validTopo);
        D_vt = D(validTopo); C_vt = C(validTopo);

        Dstar_vt = (D_vt - min(D_vt)) / max(range(D_vt), eps);
        Cstar_vt = (C_vt - min(C_vt)) / max(range(C_vt), eps);

        J_hier = 0.5*Dstar_vt + 0.5*Cstar_vt;
        [Jmin_hier, idxHier] = min(J_hier);
        alphaOptimal_hier = alphaValidTopo(idxHier);

        fprintf('Valid alphas (T=0): %s\n', mat2str(alphaValidTopo));
        fprintf('alpha* (hierarchical) = %.2f (J=%.3f)\n', alphaOptimal_hier, Jmin_hier);

        allResults.(fieldName).alphaOptimal_hier = alphaOptimal_hier;

        % --- Sensitivity to the relative weight of D and C ---
        wD_grid2D = 0.1:0.05:0.9;
        alphaStarList_hier = nan(size(wD_grid2D));

        for wi = 1:numel(wD_grid2D)
            wDi = wD_grid2D(wi);
            wCi = 1 - wDi;
            Jw = wDi*Dstar_vt + wCi*Cstar_vt;
            [~, idxW] = min(Jw);
            alphaStarList_hier(wi) = alphaValidTopo(idxW);
        end

        modeAlphaHier = mode(alphaStarList_hier);
        fracStableHier = mean(alphaStarList_hier == modeAlphaHier);
        rangeAlphaHier = [min(alphaStarList_hier), max(alphaStarList_hier)];

        fprintf('\nHierarchical criterion, weight sensitivity (wD from 0.1 to 0.9, %d combinations):\n', numel(wD_grid2D));
        fprintf('  most frequent alpha*: %.2f (in %.1f%% of the combinations)\n', modeAlphaHier, 100*fracStableHier);
        fprintf('  observed range of alpha*: [%.2f, %.2f]\n', rangeAlphaHier(1), rangeAlphaHier(2));

        allResults.(fieldName).fracStableHier = fracStableHier;
        allResults.(fieldName).rangeAlphaHier = rangeAlphaHier;

        figure(findobj('Name', 'D/C sensitivity - hierarchical criterion'));
        subplot(1, numel(sites), s);
        plot(wD_grid2D, alphaStarList_hier, 'o-', 'LineWidth', 1.5);
        xlabel('Weight w_D (fidelity), w_C = 1-w_D');
        ylabel('alpha* (hierarchical)');
        title(sprintf('%s: D/C sensitivity (%.0f%% at the mode)', siteName, 100*fracStableHier));
        grid on;
    end
end

%% ============================================================
%  LOCAL FUNCTIONS
%  ============================================================
function D = pairwiseDistances_notoolbox(A, B)
    D = sqrt(sum((permute(A,[1 3 2]) - permute(B,[3 1 2])).^2, 3));
end

function [components, sizes] = conncomp3d_notoolbox(mask)
    [nx, ny, nz] = size(mask);
    visited = false(size(mask));
    idxList = find(mask);
    components = {}; sizes = [];
    for ii = 1:numel(idxList)
        idx = idxList(ii);
        if visited(idx), continue; end
        queue = idx; visited(idx) = true; comp = idx; qi = 1;
        while qi <= numel(queue)
            cur = queue(qi); qi = qi + 1;
            [i,j,k] = ind2sub(size(mask), cur);
            for di=-1:1
                for dj=-1:1
                    for dk=-1:1
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
