%% 3_3_2_micio_extra_holdouts.m
% Additional held-out sectors for Micio Conti: slabs along the axis of the tube.
% For each slab: held-out = TLS points in the slab; constraint = all the other points;
% evaluation rectangle = the slab (no constraint point inside the slab).
% This partition differs from the main one (constraint share 70-80% instead of 33%).
%
% Two modes (variable 'modo'):
%   'lista'  -> without inversions, for each candidate slab computes the held-out points,
%               the constraint share, the evaluation elements (exclusion 1.5 m from the
%               constraint points, within 3 m of the TLS) and the cavity elements. Also
%               prints the main axis of the tube (PCA) to choose 'asse'. Row 1 is the
%               current partition (check: 355 elements, 68 in the cavity).
%   'scrivi' -> for each slab in 'scelte' creates a folder with grotta_constraint_only.txt,
%               grotta_heldout.txt, holdout_region.mat and holdout_region_ok.mat.
%               holdout_region.mat is read by 3_2_4_2_cavity_to_rect_anomalies.m, to be
%               run inside the folder of the slab.
%
% Selection rule (set before any inversion): slabs of fixed width, not overlapping the
% current partition nor each other, with at least 'minCavElems' cavity elements.

clear; clc;

%% ---------- SETTINGS ----------
modo  = 'scrivi';                        % 'lista' or 'scrivi'
asse  = 'X';                             % 'X' or 'Y' (main axis printed by the PCA)
sliceWidth = 12;                         % m
stepS = 3;                               % m, window step
minCavElems = 28;                        % minimum number of cavity elements
tlsFullFile    = 'Micio Conti_cropped_ERTarea.txt';
constraintFile = 'Micio Conti_constraint_only.txt';   % current partition (only for the check row)
holdoutFile    = 'Micio Conti_holdout_region.mat';     % current partition (only for the check row)
modelFile_U    = 'Model_Micio Conti_U.txt';            % model file (only the element coordinates are used)
alphaRadius = 1.35;
polyXY = [510340.5 4157110; 510311.1 4157106; 510290.7 4157124; 510325.0 4157133];
polyERT = polyXY;
AOI.Zmin = 251.50; AOI.Zmax = 276.93;
bufferRadius = 1.5; maxDistTLS = 3;

%% ---------- DATA ----------
tls = readmatrix(tlsFullFile); tls = tls(all(isfinite(tls(:,1:3)),2),1:3);
d = readmatrix(modelFile_U); XYZ = d(:,1:3);
inA = inpolygon(XYZ(:,1),XYZ(:,2),polyXY(:,1),polyXY(:,2)) & XYZ(:,3)>=AOI.Zmin & XYZ(:,3)<=AOI.Zmax;
XYZ = XYZ(inA,:);
shp = alphaShape(tls(:,1),tls(:,2),tls(:,3),alphaRadius);
inCavAll = inShape(shp, XYZ(:,1), XYZ(:,2), XYZ(:,3));
pad = 5;                                 % m, margin on the other axis (the slab spans the whole width)
B = [min(polyXY(:,1))-pad max(polyXY(:,1))+pad min(polyXY(:,2))-pad max(polyXY(:,2))+pad];
if strcmp(asse,'X'), col = 1; else, col = 2; end

% Slabs to write: 12 m windows with ptsOutsideZ = 0, not overlapping each other.
% Offsets (m) from the minimum X of the TLS cloud.
lo0s = min(tls(:,col));
scelte = lo0s + [18 30; 30 42];

% main axis of the tube
c = cov(tls(:,1:2)); [V,Dv] = eig(c); [~,im] = max(diag(Dv)); v = V(:,im);
fprintf('TLS points: %d ; elements inside the AOI: %d\n', size(tls,1), size(XYZ,1));
fprintf('Main axis of the tube (XY): direction (%.2f, %.2f) -> angle %.0f deg from X. TLS extent: X %.1f-%.1f, Y %.1f-%.1f\n', ...
    v(1), v(2), atan2d(v(2),v(1)), min(tls(:,1)),max(tls(:,1)),min(tls(:,2)),max(tls(:,2)));

mk = @(lo,hi) makeBox(col,lo,hi,B,AOI);
stat = @(H) sliceStats(H,tls,XYZ,inCavAll,AOI,bufferRadius,maxDistTLS);

switch modo
    case 'lista'
        % check row: current partition (constraint from file)
        L = load(holdoutFile,'HELDOUT'); cpA = readmatrix(constraintFile); cpA = cpA(:,1:3);
        sA = sliceStats(L.HELDOUT,tls,XYZ,inCavAll,AOI,bufferRadius,maxDistTLS,cpA);
        rows = [NaN NaN sA];
        lo0 = min(tls(:,col)); hi0 = max(tls(:,col));
        for a = lo0:stepS:(hi0-sliceWidth)
            rows(end+1,:) = [a a+sliceWidth stat(mk(a,a+sliceWidth))]; %#ok<SAGROW>
        end
        T = array2table(rows,'VariableNames',{'min','max','ptsHeldout','ptsConstraint','pctConstraint','ptsOutsideZ','nEval','nCavity','pctCavity'});
        writetable(T,'candidate_slices_Micio.csv'); disp(T);
        fprintf('\nRow 1 = current partition (check: 355 elements, 68 in the cavity).\n');
        ok = T(2:end,:); ok = ok(ok.nCavity>=minCavElems,:);
        fprintf('Slabs with >= %d cavity elements: %d\n', minCavElems, height(ok));
        disp(sortrows(ok,'nCavity','descend'));
        fprintf('\nSlabs in ''scelte'':\n'); R2=[];
        for i=1:size(scelte,1), R2(end+1,:)=[scelte(i,:), stat(mk(scelte(i,1),scelte(i,2)))]; end %#ok<SAGROW>
        T2=array2table(R2,'VariableNames',T.Properties.VariableNames); disp(T2); format long g; disp(scelte); format short

    case 'scrivi'
        for i = 1:size(scelte,1)
            H = mk(scelte(i,1),scelte(i,2));
            tag = sprintf('slice_%s%d_%d', asse, round(scelte(i,1)), round(scelte(i,2)));
            if ~exist(tag,'dir'), mkdir(tag); end
            isHO = inBox(tls,H);
            writematrix(tls(~isHO,:),fullfile(tag,'grotta_constraint_only.txt'),'Delimiter',' ');
            writematrix(tls(isHO,:), fullfile(tag,'grotta_heldout.txt'),'Delimiter',' ');
            HELDOUT = H;
            % holdout_region.mat: polyXY and AOI describe the HELDOUT region (rectangle + Z)
            polyXY = [H.Xmin H.Ymin; H.Xmax H.Ymin; H.Xmax H.Ymax; H.Xmin H.Ymax]; AOI.Zmin = H.Zmin; AOI.Zmax = H.Zmax; %#ok<NASGU>
            save(fullfile(tag,'holdout_region.mat'),'HELDOUT','polyXY','AOI'); save(fullfile(tag,'holdout_region_ok.mat'),'HELDOUT','polyXY','AOI');
            polyXY = polyERT;
            s = stat(H);
            fprintf('%s: held-out %d points, constraint %d (%.1f%%), evaluation elements %d (%d in the cavity)\n', tag, s(1), s(2), s(3), s(5), s(6));
            if s(4) > 0, warning('%s: %d points of the slab are outside the Z range: they remain in the constraint.', tag, s(4)); end
        end
        fprintf('\nNext, in each folder: 3_2_4_2_cavity_to_rect_anomalies.m (alphaRadius 1.35, voxelSize 0.5, rho 75000)\n');
        fprintf('-> ERTLab boxes -> H inversion -> evaluation.\n');
end

%% ---------- LOCAL FUNCTIONS ----------
function H = makeBox(col,lo,hi,B,AOI)
    H.Xmin = B(1); H.Xmax = B(2); H.Ymin = B(3); H.Ymax = B(4);
    if col == 1, H.Xmin = lo; H.Xmax = hi; else, H.Ymin = lo; H.Ymax = hi; end
    H.Zmin = AOI.Zmin; H.Zmax = AOI.Zmax;
end
function m = inBox(P,H)
    m = P(:,1)>=H.Xmin & P(:,1)<=H.Xmax & P(:,2)>=H.Ymin & P(:,2)<=H.Ymax & P(:,3)>=H.Zmin & P(:,3)<=H.Zmax;
end
function out = sliceStats(H,tls,XYZ,inCavAll,AOI,bufferRadius,maxDistTLS,cpGiven)
    isHO = inBox(tls,H);
    HXY = H; HXY.Zmin = -Inf; HXY.Zmax = Inf;
    ptsOutsideZ = nnz(inBox(tls,HXY) & ~isHO);
    if nargin >= 8, cp = cpGiven; else, cp = tls(~isHO,:); end
    rect = inBox(XYZ,H); P = XYZ(rect,:);
    if isempty(cp) || isempty(P), out = [nnz(isHO) size(cp,1) 100*size(cp,1)/size(tls,1) ptsOutsideZ 0 0 0]; return; end
    [~,dC] = knnsearch(cp,P,'K',1); [~,dT] = knnsearch(tls,P,'K',1);
    ev = (dC > bufferRadius) & (dT <= maxDistTLS);
    cav = inCavAll(rect); cav = cav(ev);
    out = [nnz(isHO), size(cp,1), 100*size(cp,1)/size(tls,1), ptsOutsideZ, nnz(ev), nnz(cav), 100*nnz(cav)/max(nnz(ev),1)];
end
