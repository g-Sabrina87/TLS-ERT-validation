%% 8_5_constraint_mesh_discretization_check.m
% Approximate estimate of how much of the constraint volume (fine-resolution boxes) is
% retained when assigned to larger cells such as those of the ERTLab mesh. The real
% ERTLab mesh is tetrahedral and non-uniform: a regular grid of side "meshSize" is an
% approximation that gives a range, not a reconstruction of the actual mesh.

clear; clc; rng(1);

%% ---------- PARAMETERS ----------
boxFile  = 'Micio_Conti_cavity_anomalies_ERTLab_75000_C.txt';  % constraint box file (7 columns)
meshSize = 1.471;                          % cell side of the ERTLab Interior Grid (m)
f = 0.25;                                  % fine step used to represent the boxes (m)
nOffsets = 20;                             % random positions of the coarse grid

%% ---------- READ BOXES ----------
B = readmatrix(boxFile);                   % Xmin Xmax Ymin Ymax Zmin Zmax Rho
nB = size(B,1);
V0 = sum((B(:,2)-B(:,1)) .* (B(:,4)-B(:,3)) .* (B(:,6)-B(:,5)));
fprintf('Boxes read: %d   Total box volume: %.1f m^3\n', nB, V0);

%% ---------- FINE MASK (centres of the fine voxels inside the union of the boxes) ----------
X0 = min(B(:,1)); Y0 = min(B(:,3)); Z0 = min(B(:,5));
nx = ceil((max(B(:,2)) - X0)/f); ny = ceil((max(B(:,4)) - Y0)/f); nz = ceil((max(B(:,6)) - Z0)/f);
M = false(nx, ny, nz);
for i = 1:nB
    i1 = max(1,  floor((B(i,1) - X0)/f + 0.5) + 1); i2 = min(nx, floor((B(i,2) - X0)/f - 0.5) + 1);
    j1 = max(1,  floor((B(i,3) - Y0)/f + 0.5) + 1); j2 = min(ny, floor((B(i,4) - Y0)/f - 0.5) + 1);
    k1 = max(1,  floor((B(i,5) - Z0)/f + 0.5) + 1); k2 = min(nz, floor((B(i,6) - Z0)/f - 0.5) + 1);
    if i2 >= i1 && j2 >= j1 && k2 >= k1
        M(i1:i2, j1:j2, k1:k2) = true;
    end
end
Vfine = nnz(M) * f^3;
fprintf('Fine-mask volume: %.1f m^3 (difference from the box volume: %.1f%%)\n', ...
    Vfine, 100*(Vfine - V0)/V0);
if abs(Vfine - V0)/V0 > 0.05
    warning('The fine mask differs by more than 5%% from the box volume: the boxes may overlap or not be aligned to a common grid.');
end

%% ---------- ASSIGNMENT TO COARSE CELLS, 3 RULES ----------
idx = find(M);
[ii, jj, kk] = ind2sub(size(M), idx);
xc = X0 + (ii - 0.5)*f; yc = Y0 + (jj - 0.5)*f; zc = Z0 + (kk - 0.5)*f;

Vany = nan(nOffsets,1); Vhalf = nan(nOffsets,1); Vfull = nan(nOffsets,1);
for t = 1:nOffsets
    o = rand(1,3) * meshSize;
    ci = floor((xc - X0 - o(1))/meshSize);
    cj = floor((yc - Y0 - o(2))/meshSize);
    ck = floor((zc - Z0 - o(3))/meshSize);
    [~, ~, g] = unique([ci cj ck], 'rows');
    cnt = accumarray(g, 1);
    cover = cnt * f^3 / meshSize^3;        % fraction of the cell covered by the constraint
    Vany(t)  = numel(cover) * meshSize^3;                       % assigned if partly inside
    Vhalf(t) = nnz(cover >= 0.5)  * meshSize^3;                 % assigned if at least half inside
    Vfull(t) = nnz(cover >= 0.97) * meshSize^3;                 % assigned only if (almost) entirely inside
end

fprintf('\n===== VOLUME RETAINED on cells of %.3f m (ratio to the box volume) =====\n', meshSize);
fprintf('Rule "partly inside"   : mean %.2f  (min %.2f, max %.2f)\n', mean(Vany/V0),  min(Vany/V0),  max(Vany/V0));
fprintf('Rule "at least half"   : mean %.2f  (min %.2f, max %.2f)\n', mean(Vhalf/V0), min(Vhalf/V0), max(Vhalf/V0));
fprintf('Rule "entirely inside" : mean %.2f  (min %.2f, max %.2f)\n', mean(Vfull/V0), min(Vfull/V0), max(Vfull/V0));
