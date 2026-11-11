%% 3_2_4_2_cavity_to_rect_anomalies.m
% Converts a TLS point cloud of a cavity into a list of axis-aligned box
% anomalies in the ERTLab format:
%
%     Xmin  Xmax  Ymin  Ymax  Zmin  Zmax  resistivity (ohm-m)
%
% The volume is clipped to the area investigated by the ERT (polygon + Z range).
% Used for the geometries of C and H (step 3.2) and for the null constraints (step 4.2:
% nullMode = true, which processes all the translated clouds of 4_1 in a loop).
% If useHoldout is true, the boxes whose centre lies inside the held-out region
% (holdoutFile) are removed.
%
% Reproducibility: the box files used in the inversions of the paper are provided in the
% repository. This script regenerates the geometry with the parameters of the paper; the result
% is close to the provided files, and identical for H at Micio Conti, but not identical for the
% others (the cleaning of the cloud and the merge of the voxels into boxes may differ).

clear; clc; close all;

%% ---------- 1) PARAMETERS ----------
site = 'Micio Conti';       % 'Micio Conti' or 'Chiusazza'
useHoldout = true;          % true: H (remove the boxes inside the held-out region); false: C (full geometry)
% Cleaning of the cloud before the alpha shape.
%   H: outlier filter + decimation (both sites).
%   C, Micio Conti: raw cloud (neither step); the alpha-shape volume is then 493.2 m^3, as in the paper.
%   C, Chiusazza: decimation only, grid from the corner of the ERT area (set below).
useOutlierFilter = useHoldout;
useDecimation    = useHoldout;
% Origin of the voxel grid: false = bounding box of the cloud clipped to the ERT area (H);
% true = corner of the ERT area (polygon minimum, AOI.Zmin), as in the original C of Chiusazza.
gridFromAOI = false;
% Null constraints (step 4.2, Micio Conti): true = process the translated clouds
% TLS_NULL_NN.txt produced by 4_1, clipped to the extent of the model mesh (the nulls lie
% partly outside the ERT quadrilateral), with the same settings as H, and write
% cavity_anomalies_NULL_NN_75000.txt (input of 4_3_check_null_boxes.m).
nullMode = false;
nullDir  = 'NULL_TLS_constraints_MicioConti';
nNull    = 30;
meshDomain = [510291.5 510397.4  4157088.0 4157141.0  251.5 276.93];   % X Y Z of Model_Micio_Conti_U

switch site
    case 'Micio Conti'
        inputFile   = 'Micio Conti_constraint_only.txt';   % cloud used for the constraint
        holdoutFile = 'Micio Conti_holdout_region.mat';
        alphaRadius = 1.35;     % alpha-shape critical radius (m)
        voxelSize   = 0.5;      % voxel size (m)
        rhoCavity   = 75000;    % resistivity assigned to the cavity (ohm-m)
        % ERT area (polygon + Z range), same as in 1_crop_tls_to_ert_area.m
        polyXY = [ 510340.5  4157110;   % P1
                   510311.1  4157106;   % P2
                   510290.7  4157124;   % P3
                   510325.0  4157133 ]; % P4
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93;
        outFile = 'Micio_Conti_cavity_anomalies_ERTLab_75000_H.txt';
    case 'Chiusazza'
        inputFile   = 'Chiusazza_constraint_only.txt';
        holdoutFile = 'Chiusazza_holdout_region.mat';
        alphaRadius = 2.50;
        voxelSize   = 1.0;
        rhoCavity   = 20000;
        polyXY = [ 514191.91  4097666.70;
                   514306.99  4097666.70;
                   514306.99  4097785.75;
                   514191.91  4097785.75 ];
        AOI.Zmin = 79.82;
        AOI.Zmax = 115.71;
        outFile = 'Chiusazza_cavity_anomalies_ERTLab_20000_H.txt';
end
if ~useHoldout                      % C: boxes from the full cloud of the site
    if strcmp(site,'Micio Conti')
        useOutlierFilter = false;  useDecimation = false;   % raw cloud
        AOI.Zmax = 278.93;                                  % top of the cropped cloud
    else
        useOutlierFilter = true;   useDecimation = true;
    end
    if strcmp(site,'Chiusazza'),   gridFromAOI = true;  end % grid from the corner of the ERT area
    inputFile = [site '_cropped_ERTarea.txt'];
    outFile   = sprintf('%s_cavity_anomalies_ERTLab_%d_C.txt', strrep(site,' ','_'), rhoCavity);
end

if nullMode
    if ~strcmp(site,'Micio Conti'), error('nullMode is set up for Micio Conti.'); end
    useHoldout = false;  useOutlierFilter = true;  useDecimation = true;  gridFromAOI = false;
    polyXY = [meshDomain(1) meshDomain(3); meshDomain(2) meshDomain(3); ...
              meshDomain(2) meshDomain(4); meshDomain(1) meshDomain(4)];
    AOI.Zmin = meshDomain(5);  AOI.Zmax = meshDomain(6);
    jobs = 1:nNull;
else
    jobs = 0;
end

for iJob = jobs
if nullMode
    inputFile = fullfile(nullDir, sprintf('TLS_NULL_%02d.txt', iJob));
    outFile   = sprintf('cavity_anomalies_NULL_%02d_%d.txt', iJob, rhoCavity);
    fprintf('\n===== NULL %02d =====\n', iJob);
end

%% ---------- 2) READ POINT CLOUD ----------
data = readmatrix(inputFile);
xyz  = data(:,1:3);
fprintf('Points read: %d\n', size(xyz,1));

%% ---------- 3) CLEANING / DECIMATION ----------
% (a) Statistical outlier removal (requires knnsearch; skipped if unavailable)
if useOutlierFilter
  try
    kNN = 8;                 % number of neighbours
    stdFactor = 2.5;         % threshold: mean + stdFactor*std of the mean neighbour distance
    [~, D] = knnsearch(xyz, xyz, 'K', kNN+1);   % column 1 is the point itself
    meanDist = mean(D(:,2:end), 2);
    thr = mean(meanDist) + stdFactor*std(meanDist);
    keep = meanDist <= thr;
    fprintf('Outlier filter: removed %d points out of %d\n', sum(~keep), numel(keep));
    xyz = xyz(keep, :);
  catch
    warning('knnsearch not available: outlier filter skipped.');
  end
else
    fprintf('Outlier filter not applied.\n');
end

% (b) Grid downsampling (mean of the points in each cell)
if useDecimation
    cellSize = voxelSize/2;
    ijk = round(xyz ./ cellSize);
    [~, ~, ic] = unique(ijk, 'rows');
    xC = accumarray(ic, xyz(:,1), [], @mean);
    yC = accumarray(ic, xyz(:,2), [], @mean);
    zC = accumarray(ic, xyz(:,3), [], @mean);
    xyz = [xC, yC, zC];
else
    fprintf('Decimation not applied.\n');
end

fprintf('Points after cleaning/decimation: %d\n', size(xyz,1));

%% ---------- 4) ALPHA SHAPE ----------
shp = alphaShape(xyz(:,1), xyz(:,2), xyz(:,3), alphaRadius);

vol = volume(shp);
fprintf('Volume enclosed by the alpha shape: %.3f m^3\n', vol);
if vol < 1e-6
    warning('The enclosed volume is almost zero: increase alphaRadius.');
end

%% ---------- 5) UNIFORM GRID AND MULTI-AXIS MERGE (limited to the ERT area) ----------
capBounds = [min(xyz(:,1)) max(xyz(:,1)) ...
             min(xyz(:,2)) max(xyz(:,2)) ...
             min(xyz(:,3)) max(xyz(:,3))];

polyBounds = [min(polyXY(:,1)) max(polyXY(:,1)) ...
              min(polyXY(:,2)) max(polyXY(:,2))];

box0 = [max(polyBounds(1), capBounds(1)), min(polyBounds(2), capBounds(2)), ...
        max(polyBounds(3), capBounds(3)), min(polyBounds(4), capBounds(4)), ...
        max(AOI.Zmin, capBounds(5)), min(AOI.Zmax, capBounds(6))];

if any(box0([2 4 6]) <= box0([1 3 5]))
    error('The intersection between the cavity and the ERT area is empty.');
end

if gridFromAOI
    box0 = [polyBounds(1) polyBounds(2) polyBounds(3) polyBounds(4) AOI.Zmin AOI.Zmax];
end
fprintf('Grid domain (bounding box of the polygon, cavity clipped):\n');
fprintf('  X: [%.2f, %.2f]  Y: [%.2f, %.2f]  Z: [%.2f, %.2f]\n', box0);

xv = box0(1) : voxelSize : box0(2);
yv = box0(3) : voxelSize : box0(4);
zv = box0(5) : voxelSize : box0(6);

[XX, YY, ZZ] = ndgrid(xv, yv, zv);
fprintf('Candidate voxels in the bounding box: %d\n', numel(XX));

% Occupancy: inside the cavity (inShape) AND inside the ERT polygon (inpolygon)
insideCavity = inShape(shp, XX(:), YY(:), ZZ(:));
insidePoly   = inpolygon(XX(:), YY(:), polyXY(:,1), polyXY(:,2));
inside = insideCavity & insidePoly;
inside = reshape(inside, size(XX));
fprintf('Occupied voxels (inside the cavity and the ERT polygon): %d\n', nnz(inside));

if nnz(inside) == 0
    warning('No voxel of the cavity falls inside the ERT area.');
end

boxes = merge_voxels_to_boxes(inside, xv, yv, zv, voxelSize);

fprintf('Number of boxes (after X/Y/Z merge): %d\n', size(boxes,1));

%% ---------- 6) REMOVAL OF THE BOXES INSIDE THE HELD-OUT REGION ----------
if useHoldout && exist(holdoutFile,'file')

    loaded = load(holdoutFile);
    disp(fieldnames(loaded)');                     % fields found in the .mat file

    % box centres
    xC = (boxes(:,1)+boxes(:,2))/2;
    yC = (boxes(:,3)+boxes(:,4))/2;
    zC = (boxes(:,5)+boxes(:,6))/2;

    if isfield(loaded,'polyXY') && isfield(loaded,'AOI')       % polygon + Z range
        insideXY = inpolygon(xC, yC, loaded.polyXY(:,1), loaded.polyXY(:,2));
        insideZ  = zC >= loaded.AOI.Zmin & zC <= loaded.AOI.Zmax;
    elseif isfield(loaded,'HELDOUT')                          % rectangle
        H = loaded.HELDOUT;
        insideXY = xC >= H.Xmin & xC <= H.Xmax & yC >= H.Ymin & yC <= H.Ymax;
        insideZ  = zC >= H.Zmin & zC <= H.Zmax;
    else
        error('%s contains neither polyXY/AOI nor HELDOUT.', holdoutFile);
    end
    overlapMask = insideXY & insideZ;

    fprintf('Held-out mode: removing %d boxes out of %d inside the excluded region.\n', ...
        sum(overlapMask), size(boxes,1));

    boxes = boxes(~overlapMask,:);
    fprintf('Remaining boxes: %d\n', size(boxes,1));

else
    fprintf('Held-out filter not applied (useHoldout = false or %s not found).\n', holdoutFile);
end

%% ---------- 7) WRITE THE ERTLab FILE ----------
fid = fopen(outFile, 'w');
for i = 1:size(boxes,1)
    fprintf(fid, '%.3f %.3f %.3f %.3f %.3f %.3f %.1f\n', ...
        boxes(i,1), boxes(i,2), boxes(i,3), ...
        boxes(i,4), boxes(i,5), boxes(i,6), rhoCavity);
end
fclose(fid);

fprintf('File written: %s\n', outFile);

% Total volume of the boxes
B = readmatrix(outFile);
V = (B(:,2)-B(:,1)).*(B(:,4)-B(:,3)).*(B(:,6)-B(:,5));
fprintf('Boxes: %d ; total volume = %.1f m^3 ; rho = %g\n', size(B,1), sum(V), unique(B(:,7)));

end   % loop over the nulls (a single pass if nullMode = false)

%% ============================================================
function boxes = merge_voxels_to_boxes(occ, xv, yv, zv, vs)
% Merges adjacent occupied voxels along X, Y and Z (several passes) into
% axis-aligned boxes.
% occ: logical array [nx x ny x nz]; xv, yv, zv: coordinates of the lower
% corner of each voxel; vs: voxel size (m)

    idx = find(occ);
    [I, J, K] = ind2sub(size(occ), idx);
    n = numel(I);

    boxes = zeros(n, 6);
    boxes(:,1) = xv(I)';        boxes(:,2) = xv(I)' + vs;
    boxes(:,3) = yv(J)';        boxes(:,4) = yv(J)' + vs;
    boxes(:,5) = zv(K)';        boxes(:,6) = zv(K)' + vs;

    % Merge along X, then Y, then Z, repeated twice
    for pass = 1:2
        boxes = merge_along_dim(boxes, 1, 2); % X
        boxes = merge_along_dim(boxes, 3, 4); % Y
        boxes = merge_along_dim(boxes, 5, 6); % Z
    end
end

function merged = merge_along_dim(boxes, loCol, hiCol)
% Merges adjacent boxes along the dimension [loCol hiCol], when the extent in
% the other 4 columns is the same (numerical tolerance tol).
    if isempty(boxes)
        merged = boxes;
        return;
    end

    tol = 1e-6;
    otherCols = setdiff(1:6, [loCol hiCol]);

    sortKey = [boxes(:, otherCols), boxes(:, loCol)];
    [~, sortIdx] = sortrows(sortKey);
    boxesSorted = boxes(sortIdx, :);

    n = size(boxesSorted, 1);
    merged = zeros(0, 6);
    i = 1;
    while i <= n
        curr = boxesSorted(i, :);
        j = i + 1;
        while j <= n && ...
              all(abs(boxesSorted(j, otherCols) - curr(otherCols)) < tol) && ...
              abs(boxesSorted(j, loCol) - curr(hiCol)) < tol
            curr(hiCol) = boxesSorted(j, hiCol);
            j = j + 1;
        end
        merged(end+1, :) = curr; %#ok<AGROW>
        i = j;
    end
end
