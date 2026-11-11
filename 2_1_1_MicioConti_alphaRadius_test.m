%% 2_1_1_MicioConti_alphaRadius_test.m
% Elbow method and geometric index J(alpha) for the choice of the alpha-shape
% radius at Micio Conti. J(alpha) = fidelity (bidirectional Chamfer distance)
% + topological fragmentation penalty (alpha-shape regions)
% + voxel-coverage fragmentation penalty (voxel size = size used for the boxes).

clear; clc; close all;

%% ---------- 1) PARAMETERS ----------
tlsFile = 'Micio Conti_cropped_ERTarea.txt';   % TLS point cloud of Micio Conti

alphaTest = 0.75:0.1:2.0;   % alphaRadius values to test (m)

% ERT area (polygon)
polyXY = [ 510340.5  4157110;   % P1
           510311.1  4157106;   % P2
           510290.7  4157124;   % P3
           510325.0  4157133 ]; % P4
AOI.Zmin = 251.50;
AOI.Zmax = 276.93;

voxelSize = 0.5;   % same value used in 3_2_4_2_cavity_to_rect_anomalies.m

lambdaFragmentation = 0.5;   % weight of the topological penalty (alpha-shape regions)
lambdaCoverage      = 2.0;   % weight of the voxel-coverage fragmentation penalty

%% ---------- 2) READ POINT CLOUD ----------
data = readmatrix(tlsFile);
xyz  = data(:,1:3);
fprintf('Points read: %d\n', size(xyz,1));

%% ---------- 3) GRID DOMAIN (bounding box of the polygon, as in 3_2_4_2_cavity_to_rect_anomalies.m) ----------
capBounds  = [min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) min(xyz(:,3)) max(xyz(:,3))];
polyBounds = [min(polyXY(:,1)) max(polyXY(:,1)) min(polyXY(:,2)) max(polyXY(:,2))];

box0 = [max(polyBounds(1), capBounds(1)), min(polyBounds(2), capBounds(2)), ...
        max(polyBounds(3), capBounds(3)), min(polyBounds(4), capBounds(4)), ...
        max(AOI.Zmin, capBounds(5)), min(AOI.Zmax, capBounds(6))];

xv = box0(1):voxelSize:box0(2);
yv = box0(3):voxelSize:box0(4);
zv = box0(5):voxelSize:box0(6);
[XX, YY, ZZ] = ndgrid(xv, yv, zv);
insidePolyMask = inpolygon(XX(:), YY(:), polyXY(:,1), polyXY(:,2));
insidePolyMask = reshape(insidePolyMask, size(XX));

fprintf('Grid generated: %d candidate voxels inside the ERT area (polygon)\n', nnz(insidePolyMask));

%% ---------- 4) COMPUTATION FOR EACH ALPHARADIUS ----------
nTest = numel(alphaTest);
volTest = nan(nTest,1);
nRegions = nan(nTest,1);
fidelityDist = nan(nTest,1);
coverageFragFraction = nan(nTest,1);
Jvalues = nan(nTest,1);

maxPointsForFidelity = 5000;
if size(xyz,1) > maxPointsForFidelity
    rng(1);
    idxSample = randperm(size(xyz,1), maxPointsForFidelity)';
    xyzForFidelity = xyz(idxSample,:);
else
    xyzForFidelity = xyz;
end

fprintf('\nComputing J(alpha) for %d values...\n', nTest);
for i = 1:nTest
    shpTest = alphaShape(xyz(:,1), xyz(:,2), xyz(:,3), alphaTest(i));
    volTest(i) = volume(shpTest);
    try
        nRegions(i) = numRegions(shpTest);
    catch
        nRegions(i) = NaN;
    end

    if volTest(i) < 1e-6
        fidelityDist(i) = NaN; coverageFragFraction(i) = NaN; Jvalues(i) = Inf;
        fprintf('  alpha=%.2f -> degenerate shape, discarded\n', alphaTest(i));
        continue;
    end

    % --- Term 1: fidelity (bidirectional Chamfer distance) ---
    try
        [tri, triPts] = boundaryFacets(shpTest);
    catch
        tri = [];
    end
    if isempty(tri)
        fidelityDist(i) = NaN; coverageFragFraction(i) = NaN; Jvalues(i) = Inf;
        continue;
    end
    nSamplesPerTri = 3;
    v1 = triPts(tri(:,1),:); v2 = triPts(tri(:,2),:); v3 = triPts(tri(:,3),:);
    surfSamples = zeros(size(tri,1)*nSamplesPerTri, 3);
    k = 0;
    for t = 1:size(tri,1)
        for s = 1:nSamplesPerTri
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
    fidelityDist(i) = mean([mean(min(d1,[],2)), mean(min(d1,[],1))]);

    % --- Term 2: topological fragmentation of the alpha shape ---
    fragPenaltyTopo = lambdaFragmentation * max(0, nRegions(i) - 1);

    % --- Term 3: continuity of the voxel coverage inside the ERT area ---
    insideCavity = inShape(shpTest, XX(:), YY(:), ZZ(:));
    insideCavity = reshape(insideCavity, size(XX));
    occMask = insideCavity & insidePolyMask;

    nOccupied = nnz(occMask);
    if nOccupied == 0
        coverageFragFraction(i) = 1;
    else
        try
            CC = bwconncomp(occMask, 26);
            compSizes = cellfun(@numel, CC.PixelIdxList);
        catch
            [~, compSizes] = conncomp3d_notoolbox(occMask);
        end
        mainSize = max(compSizes);
        coverageFragFraction(i) = 1 - mainSize/nOccupied;   % 0 = fully connected
    end
    fragPenaltyCoverage = lambdaCoverage * coverageFragFraction(i);

    Jvalues(i) = fidelityDist(i) + fragPenaltyTopo + fragPenaltyCoverage;

    fprintf(['  alpha=%.2f -> fidelity=%.3f m, alpha-shape regions=%d, ' ...
             'coverage fragmentation=%.1f%%, J=%.3f\n'], ...
        alphaTest(i), fidelityDist(i), nRegions(i), 100*coverageFragFraction(i), Jvalues(i));
end

[Jmin, bestIdx] = min(Jvalues);
alphaOptimal = alphaTest(bestIdx);

fprintf('\n=== ALPHARADIUS AT THE MINIMUM OF J ===\n');
fprintf('alphaRadius = %.2f (J=%.3f, fidelity=%.3f m, regions=%d, coverage fragmentation=%.1f%%, volume=%.1f m^3)\n', ...
    alphaOptimal, Jmin, fidelityDist(bestIdx), nRegions(bestIdx), 100*coverageFragFraction(bestIdx), volTest(bestIdx));

%% ---------- 5) PLOTS ----------
figure;
subplot(2,1,1);
plot(alphaTest, volTest, 'o-', 'LineWidth', 1.5);
xlabel('alphaRadius (m)'); ylabel('Volume (m^3)');
title('Elbow method - Micio Conti'); grid on;

subplot(2,1,2);
plot(alphaTest, nRegions, 's-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
xlabel('alphaRadius (m)'); ylabel('Alpha-shape regions');
title('Topological fragmentation'); grid on;

figure;
subplot(2,1,1);
plot(alphaTest, Jvalues, 'o-', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.7]);
hold on;
plot(alphaOptimal, Jmin, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
xlabel('alphaRadius (m)'); ylabel('J(alpha)');
title(sprintf('J(alpha): minimum at alpha=%.2f', alphaOptimal));
legend('J(alpha)', 'Minimum', 'Location', 'best'); grid on;

subplot(2,1,2);
yyaxis left;
plot(alphaTest, fidelityDist, 's-', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);
ylabel('Fidelity (m)');
yyaxis right;
plot(alphaTest, 100*coverageFragFraction, '^-', 'Color', [0.8 0.3 0.1], 'LineWidth', 1.5);
ylabel('Coverage fragmentation (%)');
xlabel('alphaRadius (m)');
title('Components: fidelity vs voxel-coverage fragmentation');
legend('Fidelity (Chamfer, m)', 'Coverage fragmentation (%)', 'Location', 'best');
grid on;

%% ---------- 6) VOXEL COVERAGE AT alphaOptimal ----------
shpBest = alphaShape(xyz(:,1), xyz(:,2), xyz(:,3), alphaOptimal);
insideCavityBest = inShape(shpBest, XX(:), YY(:), ZZ(:));
insideCavityBest = reshape(insideCavityBest, size(XX));
occMaskBest = insideCavityBest & insidePolyMask;

figure;
plot(shpBest, 'FaceColor', [0.9 0.9 0.6], 'FaceAlpha', 0.3, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.2);
hold on;
occIdx = find(occMaskBest);
scatter3(XX(occIdx), YY(occIdx), ZZ(occIdx), 15, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
axis equal; grid on; view(2);
title(sprintf('Voxel coverage at alphaRadius=%.2f (fragmentation=%.1f%%)', ...
    alphaOptimal, 100*coverageFragFraction(bestIdx)));
xlabel('X'); ylabel('Y');

%% ============================================================
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
