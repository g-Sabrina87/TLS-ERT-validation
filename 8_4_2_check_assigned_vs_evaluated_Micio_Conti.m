%% 8_4_2_check_assigned_vs_evaluated_Micio_Conti.m
% Check: do the elements ASSIGNED in the initial model of H (those carrying the imposed
% value) fall among the EVALUATED elements? The exclusion from the evaluation is defined
% by the constraint POINTS (radius 1.5 m), not by the assigned elements; this script
% verifies that the two sets do not overlap.
%
% The assigned elements are read from the initial model of ERTLab+ (hexahedral mesh with
% the resistivity of each element): the elements whose resistivity equals the imposed
% value. The same set is obtained by assigning to the mesh the blocks of H (rule: the
% CENTROID of the element falls inside a block); the script prints both counts.
% Block file (H): Xmin Xmax Ymin Ymax Zmin Zmax rho.

clear; clc;

boxFileH         = 'Micio_Conti_cavity_anomalies_ERTLab_75000_H.txt';   % blocks of H (constraint region only)
initialModelFile = 'initial_model-holdout_75000.txt';   % initial model of H (ERTLab+ mesh); '' to skip
rhoImposed       = 75000;
modelFile      = 'Model_Micio Conti_U.txt';        % U: used only for the coordinates of the evaluated elements
holdoutFile    = 'Micio Conti_holdout_region.mat';
constraintFile = 'Micio Conti_constraint_only.txt';
tlsFullFile    = 'Micio Conti_cropped_ERTarea.txt';
polyXY = [ 510340.5  4157110;
           510311.1  4157106;
           510290.7  4157124;
           510325.0  4157133 ];
AOI.Zmin = 251.50; AOI.Zmax = 276.93;
bufferRadius = 1.5; maxDistTLS = 3;

load(holdoutFile, 'HELDOUT');
B = readmatrix(boxFileH); B = B(all(isfinite(B(:,1:6)),2), :);
d = readmatrix(modelFile); P = d(:,1:3);
cp  = readmatrix(constraintFile); cp  = cp(all(isfinite(cp(:,1:3)),2), 1:3);
tls = readmatrix(tlsFullFile);    tls = tls(all(isfinite(tls(:,1:3)),2), 1:3);

%% ---------- ASSIGNED ELEMENTS ----------
if ~isempty(initialModelFile) && isfile(initialModelFile)
    [Cm, rhoEl] = read_ertlab_initial_model(initialModelFile);   % centroids and resistivity of the elements
    isAss = abs(rhoEl - rhoImposed) < 1e-6*rhoImposed;
    Pass  = Cm(isAss,:);
    inBox = false(size(Cm,1),1);
    for b = 1:size(B,1)
        inBox = inBox | (Cm(:,1)>=B(b,1) & Cm(:,1)<=B(b,2) & Cm(:,2)>=B(b,3) & Cm(:,2)<=B(b,4) & Cm(:,3)>=B(b,5) & Cm(:,3)<=B(b,6));
    end
    fprintf('Blocks: %d. Elements of the initial model at %g ohm-m: %d ; centroid rule on the mesh: %d (identical sets: %d)\n', ...
        size(B,1), rhoImposed, nnz(isAss), nnz(inBox), isequal(inBox, isAss));
else
    % fallback: centroid rule on the coordinates of the exported model (these have a Z about
    % 0.7 m lower than the centroids of the ERTLab+ mesh: fewer elements are found)
    warning('Initial model not found: using the exported coordinates (approximate).');
    inBox = false(size(P,1),1);
    for b = 1:size(B,1)
        inBox = inBox | (P(:,1)>=B(b,1) & P(:,1)<=B(b,2) & P(:,2)>=B(b,3) & P(:,2)<=B(b,4) & P(:,3)>=B(b,5) & P(:,3)<=B(b,6));
    end
    Pass = P(inBox,:);
    fprintf('Blocks: %d. Assigned elements (centroid in a block, exported coordinates): %d\n', size(B,1), size(Pass,1));
end

%% ---------- EVALUATION SET (same selection as 7_2_null_diagnostics_buffer.m) ----------
inArea = inpolygon(P(:,1),P(:,2),polyXY(:,1),polyXY(:,2)) & P(:,3)>=AOI.Zmin & P(:,3)<=AOI.Zmax;
inHO = P(:,1)>=HELDOUT.Xmin & P(:,1)<=HELDOUT.Xmax & P(:,2)>=HELDOUT.Ymin & P(:,2)<=HELDOUT.Ymax & ...
       P(:,3)>=HELDOUT.Zmin & P(:,3)<=HELDOUT.Zmax;
rect = inArea & inHO;
[~, dC] = knnsearch(cp,  P(rect,:), 'K', 1);
[~, dT] = knnsearch(tls, P(rect,:), 'K', 1);
evalSet = false(size(P,1),1);
evalSet(find(rect)) = (dC > bufferRadius) & (dT <= maxDistTLS); %#ok<FNDSB>
fprintf('Elements in the held-out rectangle: %d ; in the evaluation set: %d\n', nnz(rect), nnz(evalSet));

inHOass = Pass(:,1)>=HELDOUT.Xmin & Pass(:,1)<=HELDOUT.Xmax & Pass(:,2)>=HELDOUT.Ymin & Pass(:,2)<=HELDOUT.Ymax & ...
          Pass(:,3)>=HELDOUT.Zmin & Pass(:,3)<=HELDOUT.Zmax;
fprintf('\nASSIGNED elements: %d ; inside the held-out rectangle: %d\n', size(Pass,1), nnz(inHOass));

% distance of the evaluated elements from the assigned ones
[~, dA] = knnsearch(Pass, P(evalSet,:), 'K', 1);
fprintf('Distance of the evaluated elements from the nearest assigned element: min %.2f, median %.2f m\n', min(dA), median(dA));
fprintf('  evaluated elements within 1.5 m of an assigned element: %d of %d\n', nnz(dA <= 1.5), numel(dA));

%% ============================================================
function [C, rho] = read_ertlab_initial_model(fname)
% Reads an ERTLab+ initial model (hexahedral mesh): returns the centroid and the resistivity
% of each element. Sections: "# NN" (number of nodes), "# NE" (number of elements),
% "# Nodes coordinates", "# Elements-resistivity" (8 node indices + resistivity per element).
    fid = fopen(fname, 'r'); k = 0; hN = 0; hE = 0; NN = 0; NE = 0;
    while true
        l = fgetl(fid); if ~ischar(l), break; end
        k = k + 1; t = strtrim(l);
        if strcmp(t, '# NN'), fgetl(fid); l2 = fgetl(fid); k = k + 2; NN = str2double(l2);
        elseif strcmp(t, '# NE'), fgetl(fid); l2 = fgetl(fid); k = k + 2; NE = str2double(l2);
        elseif strncmp(t, '# Nodes coordinates', 19), hN = k;
        elseif strncmp(t, '# Elements-resistivity', 22), hE = k; break;
        end
    end
    fclose(fid);
    if NN == 0 || NE == 0 || hN == 0 || hE == 0, error('Unexpected format of %s', fname); end
    fid = fopen(fname, 'r');
    N = textscan(fid, '%f %f %f %*[^\n]', NN, 'HeaderLines', hN + 2);   % nodes: X Y Z
    fclose(fid);
    nodes = [N{1}, N{2}, N{3}];
    fid = fopen(fname, 'r');
    E = textscan(fid, '%f %f %f %f %f %f %f %f %f', NE, 'HeaderLines', hE + 1);   % 8 nodes + resistivity
    fclose(fid);
    idx = [E{1:8}]; rho = E{9};
    C = zeros(NE, 3);
    for c = 1:3
        v = nodes(:, c); C(:, c) = mean(v(idx), 2);
    end
end
