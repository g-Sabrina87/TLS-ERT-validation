%% 4_3_check_null_boxes.m
% Checks the boxes of the nulls before inverting them: number of boxes, total volume
% and whether they lie inside the mesh domain. Compared with H (reference).
%
% Box files: one row per box  Xmin Xmax Ymin Ymax Zmin Zmax rho
% (format of 3_2_4_2_cavity_to_rect_anomalies.m). Text/header rows are ignored.
% Output: check_null_boxes.csv

clear; clc;

boxDir     = '';                                        % folder of the files ('' = current folder)
boxPattern = 'cavity_anomalies_NULL_%02d_75000.txt';    % file names of the nulls
nNull      = 30;
refFile    = 'Micio_Conti_cavity_anomalies_ERTLab_75000_H.txt';   % boxes of H (reference); '' to skip
DOM = [510291.5 510397.4  4157088.0 4157141.0  251.5 276.93];  % X Y Z domain of Model_Micio_U
tol = 0.5;                                  % m, tolerance on the domain edge

rows = {};
if ~isempty(refFile) && isfile(refFile)
    [n0, v0, out0] = boxStats(refFile, DOM, tol);
    fprintf('REFERENCE H: %d boxes, total volume %.1f m^3, outside the domain: %d\n\n', n0, v0, out0);
    rows(end+1,:) = {0, n0, v0, out0}; %#ok<SAGROW>
else
    n0 = NaN; v0 = NaN;
end

for k = 1:nNull
    f = fullfile(boxDir, sprintf(boxPattern, k));
    if ~isfile(f), fprintf('NULL %02d: missing file (%s)\n', k, f); continue; end
    [n, v, nOut] = boxStats(f, DOM, tol);
    flag = '';
    if n == 0, flag = ' <-- EMPTY';
    elseif nOut > 0, flag = ' <-- BOXES OUTSIDE THE DOMAIN';
    elseif ~isnan(v0) && (v < 0.7*v0 || v > 1.3*v0), flag = ' <-- VOLUME DIFFERENT (>30%) FROM H';
    end
    fprintf('NULL %02d: %4d boxes | volume %7.1f m^3 | outside the domain: %d%s\n', k, n, v, nOut, flag);
    rows(end+1,:) = {k, n, v, nOut}; %#ok<SAGROW>
end

if isempty(rows)
    error('No file found. Current folder: %s', pwd);
end
T = cell2table(rows, 'VariableNames', {'NullID','nBoxes','volume_m3','nOutsideDomain'});
writetable(T, 'check_null_boxes.csv');
nn = T(T.NullID>0,:);
fprintf('\nNulls: boxes %d-%d (median %g) | volume %.1f-%.1f m^3 (median %.1f)\n', ...
    min(nn.nBoxes), max(nn.nBoxes), median(nn.nBoxes), min(nn.volume_m3), max(nn.volume_m3), median(nn.volume_m3));

function [n, v, nOut] = boxStats(f, DOM, tol)
    B = readmatrix(f);
    B = B(all(isfinite(B(:,1:6)),2), :);
    n = size(B,1);
    if n == 0, v = 0; nOut = 0; return; end
    v = sum((B(:,2)-B(:,1)).*(B(:,4)-B(:,3)).*(B(:,6)-B(:,5)));
    out = B(:,1) < DOM(1)-tol | B(:,2) > DOM(2)+tol | ...
          B(:,3) < DOM(3)-tol | B(:,4) > DOM(4)+tol | ...
          B(:,5) < DOM(5)-tol | B(:,6) > DOM(6)+tol;
    nOut = sum(out);
end
