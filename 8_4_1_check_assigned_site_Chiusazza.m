%% 8_4_1_check_assigned_site_Chiusazza.m
% Element-level quantification of the imposed constraint (C and H) for Chiusazza.
% (For Micio Conti see 8_4_2_check_assigned_vs_evaluated_Micio_Conti.m.)
% Assignment rule: an element is assigned if its CENTROID falls inside a box.
%
% For C and H:
%   - number of assigned elements and NOMINAL volume (n x dx*dy*dz; the mesh follows the
%     topography, so this is an approximation) compared with the box volume
%   - fraction of the box volume covered by assigned cells (nominal cells centred on the
%     centroids) and fraction of the assigned-cell volume lying outside the boxes
%   - for H: number of assigned elements falling in the evaluation set (expected: 0)
% Box files have 7 columns: Xmin Xmax Ymin Ymax Zmin Zmax rho.

clear; clc; rng(1);

%% ---------- SETTINGS (Chiusazza) ----------
boxFileC = 'Chiusazza_cavity_anomalies_ERTLab_20000_C.txt';        % boxes of C (full geometry), 20000 Ohm*m
boxFileH = 'Chiusazza_cavity_anomalies_ERTLab_20000_H.txt';                   % boxes of H (constraint region only), 20000 Ohm*m
modelFile = 'Model_Chiusazza_U.txt';                    % element coordinates only
holdoutFile = 'Chiusazza_holdout_region.mat';
constraintFile = 'Chiusazza_constraint_only.txt';
tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';
alphaRadius = 2.50;
cellDims = [1.075 1.002 1.002];                        % nominal dx dy dz of the mesh (m)
polyXY = [ 514191.91  4097666.70;
           514306.99  4097666.70;
           514306.99  4097785.75;
           514191.91  4097785.75 ];
AOI.Zmin = 79.82;
AOI.Zmax = 115.71;
bufferRadius = 1.0; maxDistTLS = 3;
f = 0.25;                                              % fine step for the box mask (m)

%% ---------- DATA ----------
d = readmatrix(modelFile); XYZ = d(:,1:3);
tls = readmatrix(tlsFullFile); tls = tls(all(isfinite(tls(:,1:3)),2),1:3);
shp = alphaShape(tls(:,1),tls(:,2),tls(:,3),alphaRadius);
fprintf('Model elements: %d ; volume of the full-TLS alpha shape (alpha %.2f): %.1f m^3\n', size(XYZ,1), alphaRadius, volume(shp));

report('C', boxFileC, XYZ, cellDims, f);
Hres = report('H', boxFileH, XYZ, cellDims, f);

%% ---------- H: ASSIGNED vs EVALUATION SET ----------
load(holdoutFile,'HELDOUT');
cp = readmatrix(constraintFile); cp = cp(all(isfinite(cp(:,1:3)),2),1:3);
inArea = inpolygon(XYZ(:,1),XYZ(:,2),polyXY(:,1),polyXY(:,2)) & XYZ(:,3)>=AOI.Zmin & XYZ(:,3)<=AOI.Zmax;
inHO = XYZ(:,1)>=HELDOUT.Xmin & XYZ(:,1)<=HELDOUT.Xmax & XYZ(:,2)>=HELDOUT.Ymin & XYZ(:,2)<=HELDOUT.Ymax & ...
       XYZ(:,3)>=HELDOUT.Zmin & XYZ(:,3)<=HELDOUT.Zmax;
rect = inArea & inHO;
[~,dC] = knnsearch(cp, XYZ(rect,:), 'K', 1); [~,dT] = knnsearch(tls, XYZ(rect,:), 'K', 1);
evalSet = false(size(XYZ,1),1); evalSet(find(rect)) = (dC > bufferRadius) & (dT <= maxDistTLS); %#ok<FNDSB>
fprintf('\nEvaluation set: %d elements (expected 11758). H elements assigned in the held-out rectangle: %d ; in the evaluation set: %d (expected 0)\n', ...
    nnz(evalSet), nnz(Hres.assigned & rect), nnz(Hres.assigned & evalSet));
if any(Hres.assigned)
    [~,dA] = knnsearch(XYZ(Hres.assigned,:), XYZ(evalSet,:), 'K', 1);
    fprintf('Distance of the evaluated elements from the nearest assigned element: min %.2f, median %.2f m\n', min(dA), median(dA));
end

%% ---------- LOCAL FUNCTION ----------
function R = report(tag, boxFile, XYZ, cd, f)
    B = readmatrix(boxFile); B = B(all(isfinite(B(:,1:6)),2),:);
    nB = size(B,1); V0 = sum((B(:,2)-B(:,1)).*(B(:,4)-B(:,3)).*(B(:,6)-B(:,5)));
    asg = false(size(XYZ,1),1);
    for b = 1:nB
        asg = asg | (XYZ(:,1)>=B(b,1) & XYZ(:,1)<=B(b,2) & XYZ(:,2)>=B(b,3) & XYZ(:,2)<=B(b,4) & XYZ(:,3)>=B(b,5) & XYZ(:,3)<=B(b,6));
    end
    n = nnz(asg); Vnom = n*prod(cd);
    % fine mask of the boxes
    X0 = min(B(:,1)); Y0 = min(B(:,3)); Z0 = min(B(:,5));
    nx = ceil((max(B(:,2))-X0)/f); ny = ceil((max(B(:,4))-Y0)/f); nz = ceil((max(B(:,6))-Z0)/f);
    M = false(nx,ny,nz);
    for b = 1:nB
        i1 = max(1,floor((B(b,1)-X0)/f+0.5)+1); i2 = min(nx,floor((B(b,2)-X0)/f-0.5)+1);
        j1 = max(1,floor((B(b,3)-Y0)/f+0.5)+1); j2 = min(ny,floor((B(b,4)-Y0)/f-0.5)+1);
        k1 = max(1,floor((B(b,5)-Z0)/f+0.5)+1); k2 = min(nz,floor((B(b,6)-Z0)/f-0.5)+1);
        if i2>=i1 && j2>=j1 && k2>=k1, M(i1:i2,j1:j2,k1:k2) = true; end
    end
    idx = find(M); [ii,jj,kk] = ind2sub(size(M),idx);
    pc = [X0+(ii-0.5)*f, Y0+(jj-0.5)*f, Z0+(kk-0.5)*f];          % centres of the fine voxels inside the boxes
    Vfine = size(pc,1)*f^3;
    % fraction of the box volume covered by assigned cells (nominal cells centred on the centroids)
    Cs = XYZ(asg,:)./cd;
    [~,dch] = knnsearch(Cs, pc./cd, 'K', 1, 'Distance', 'chebychev');
    covBox = mean(dch <= 0.5);
    % fraction of the assigned-cell volume outside the boxes (Monte Carlo sampling)
    nPer = 40; A = XYZ(asg,:); m = size(A,1);
    S = repelem(A,nPer,1) + (rand(m*nPer,3)-0.5).*cd;
    ix = floor((S(:,1)-X0)/f)+1; iy = floor((S(:,2)-Y0)/f)+1; iz = floor((S(:,3)-Z0)/f)+1;
    inb = ix>=1 & ix<=nx & iy>=1 & iy<=ny & iz>=1 & iz<=nz;
    inside = false(size(S,1),1); inside(inb) = M(sub2ind(size(M),ix(inb),iy(inb),iz(inb)));
    outFrac = 1 - mean(inside);
    fprintf('\n===== %s : %d boxes, box volume %.1f m^3 (fine mask %.1f) =====\n', tag, nB, V0, Vfine);
    fprintf('Assigned elements: %d ; nominal volume %.1f m^3 (ratio to the boxes %.2f)\n', n, Vnom, Vnom/V0);
    fprintf('Fraction of the box volume covered by assigned cells: %.2f\n', covBox);
    fprintf('Fraction of the assigned-cell volume outside the boxes: %.2f\n', outFrac);
    R.assigned = asg; R.n = n;
end
