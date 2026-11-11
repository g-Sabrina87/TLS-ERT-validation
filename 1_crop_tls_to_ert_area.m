%% 1_crop_tls_to_ert_area.m
% Crops the TLS point cloud of the cavity to the area investigated by the ERT:
% 2D polygon in plan view (XY) + Z range.

clear; clc; close all;

%% ---------- 1) PARAMETERS ----------
site = 'Micio Conti';        % 'Micio Conti' or 'Chiusazza'

switch site
    case 'Micio Conti'
        inputFile  = 'Micio Conti_TLS.txt';
        outputFile = 'Micio Conti_cropped_ERTarea.txt';
        % Polygon vertices P1 -> P2 -> P3 -> P4 (closed automatically)
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 278.93;   % crop margin above the evaluation top (276.93 m): gives the 65,742 points of the paper
    case 'Chiusazza'
        inputFile  = 'Chiusazza_TLS.txt';
        outputFile = 'Chiusazza_cropped_ERTarea.txt';
        % Exact AOI limits (same as in 2_1_2): rounding them to 2 decimals changes
        % the crop by a few points on the border (50,704 instead of the 50,705 of the paper).
        polyXY = [ 514191.9125  4097666.704;
                   514306.9915  4097666.704;
                   514306.9915  4097785.755;
                   514191.9125  4097785.755 ];
        AOI.Zmin = 79.82597;
        AOI.Zmax = 115.70548;
    otherwise
        error('Unknown site: %s', site);
end

%% ---------- 2) READ POINT CLOUD ----------
data = readmatrix(inputFile);
xyz  = data(:,1:3);
extraCols = data(:,4:end);   % any extra columns (intensity, RGB, ...)

fprintf('Points read: %d\n', size(xyz,1));

%% ---------- 3) CROP ----------
inXY = inpolygon(xyz(:,1), xyz(:,2), polyXY(:,1), polyXY(:,2));
inZ  = xyz(:,3) >= AOI.Zmin & xyz(:,3) <= AOI.Zmax;
keep = inXY & inZ;

xyzCropped   = xyz(keep,:);
extraCropped = extraCols(keep,:);

fprintf('Points inside the ERT area: %d (%.1f%% of the total)\n', ...
    sum(keep), 100*sum(keep)/numel(keep));

if sum(keep) == 0
    warning('No point falls inside the ERT area: check the TLS-ERT coregistration and the polygon.');
end

%% ---------- 4) SAVE ----------
writematrix([xyzCropped, extraCropped], outputFile, 'Delimiter', ' ');
fprintf('Cropped file written: %s\n', outputFile);
