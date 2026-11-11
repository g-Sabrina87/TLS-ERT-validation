%% 4_4_null_position_check.m
% Position of the nulls with respect to (a) the held-out sector, (b) the cavity
% (full TLS), (c) the electrodes. No inversions. A null is the translated
% constraint cloud:
%     null k = constraint cloud + [tx ty tz]_k      (columns 2-4 of the NULL_summary file)
%
% For each null (row 0 = real constraint, zero translation, as reference):
%   pctInHeldout  : % of the null points inside the held-out box (X,Y,Z)
%   pctNearCavity : % of the null points within nearRadius of the full TLS
%   minDistTLS    : minimum distance from the full TLS (at Micio Conti it must match the
%                   MinDistToTLS_m column of the file)
%   medDistElec   : median 3D distance of the null points from the nearest electrode
%   nElec10       : mean number of electrodes within 10 m (horizontal) of the null points
% Output: null_position_<site>.csv

clear; clc;

site = 'Micio Conti';      % 'Micio Conti' or 'Chiusazza'
nearRadius = 3;            % m
maxPts = 20000;            % subsampling of the constraint points (speed)

switch site
    case 'Micio Conti'
        summaryFile    = 'NULL_summary_Micio Conti.csv';              % columns: id, dx, dy, dz, ...
        constraintFile = 'Micio Conti_constraint_only.txt';
        tlsFullFile    = 'Micio Conti_cropped_ERTarea.txt';
        holdoutFile    = 'Micio Conti_holdout_region.mat';
        electrodeFile  = 'electrodes_MicioConti.csv';             % cable,id,X,Y,Z
    case 'Chiusazza'
        summaryFile    = 'NULL_summary_Chiusazza.csv';                      % id, tx, ty, tz, ...
        constraintFile = 'Chiusazza_constraint_only.txt';
        tlsFullFile    = 'Chiusazza_cropped_ERTarea.txt';
        holdoutFile    = 'Chiusazza_holdout_region.mat';
        electrodeFile  = 'electrodes_Chiusazza.csv';
end

load(holdoutFile, 'HELDOUT');
S   = readmatrix(summaryFile);    S   = S(all(isfinite(S(:,1:4)),2), :);
cp  = readmatrix(constraintFile); cp  = cp(all(isfinite(cp(:,1:3)),2), 1:3);
tls = readmatrix(tlsFullFile);    tls = tls(all(isfinite(tls(:,1:3)),2), 1:3);
El  = readmatrix(electrodeFile);  El  = El(:,3:5);
if size(cp,1) > maxPts, cp = cp(round(linspace(1,size(cp,1),maxPts)),:); end

nN = size(S,1);
R = nan(nN+1, 9);   % id tx ty tz pctInHO pctNear minDistTLS medDistElec nElec10
for k = 0:nN
    if k == 0, t = [0 0 0]; id = 0; else, t = S(k,2:4); id = S(k,1); end
    P = cp + t;
    inHO = P(:,1)>=HELDOUT.Xmin & P(:,1)<=HELDOUT.Xmax & P(:,2)>=HELDOUT.Ymin & P(:,2)<=HELDOUT.Ymax & ...
           P(:,3)>=HELDOUT.Zmin & P(:,3)<=HELDOUT.Zmax;
    [~, dT] = knnsearch(tls, P, 'K', 1);
    [~, dE] = knnsearch(El,  P, 'K', 1);
    idx = rangesearch(El(:,1:2), P(1:5:end,1:2), 10);
    R(k+1,:) = [id, t, 100*mean(inHO), 100*mean(dT <= nearRadius), min(dT), median(dE), mean(cellfun(@numel, idx))];
end

T = array2table(R, 'VariableNames', {'NullID','tx','ty','tz','pctInHeldout','pctNearCavity','minDistTLS','medDistElec','nElec10'});
writetable(T, sprintf('null_position_%s.csv', strrep(site,' ','')));
disp(T);

N = T(T.NullID>0,:);
fprintf('\n%s: %d nulls\n', site, height(N));
fprintf('  points in the held-out sector: mean %.1f%% (min %.1f, max %.1f); nulls with >10%% of the points inside: %d; with >50%%: %d\n', ...
    mean(N.pctInHeldout), min(N.pctInHeldout), max(N.pctInHeldout), sum(N.pctInHeldout>10), sum(N.pctInHeldout>50));
fprintf('  points within %g m of the cavity: mean %.1f%% (min %.1f, max %.1f); nulls with >15%%: %d\n', ...
    nearRadius, mean(N.pctNearCavity), min(N.pctNearCavity), max(N.pctNearCavity), sum(N.pctNearCavity>15));
fprintf('  median distance from the nearest electrode: nulls %.1f-%.1f m (mean %.1f) | real %.1f m\n', ...
    min(N.medDistElec), max(N.medDistElec), mean(N.medDistElec), T.medDistElec(1));
fprintf('  electrodes within 10 m (mean): nulls %.1f-%.1f (mean %.1f) | real %.1f\n', ...
    min(N.nElec10), max(N.nElec10), mean(N.nElec10), T.nElec10(1));
if strcmp(site,'Micio Conti') && size(S,2) >= 5
    fprintf('  CHECK: recomputed minDistTLS vs file column, maximum difference %.2f m\n', ...
        max(abs(N.minDistTLS - S(:,5))));
end
