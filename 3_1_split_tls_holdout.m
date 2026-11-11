%% 3_1_split_tls_holdout.m
% Splits the TLS point cloud of the cavity into two parts:
%   - CONSTRAINT: used to generate the box anomalies (ERTLab constraint)
%   - HELDOUT:    excluded from the constraint, used only to evaluate the inversion
% The held-out region is defined by an XY polygon + a Z range.

clear; clc; close all;

%% ---------- 1) PARAMETERS ----------
site = 'Micio Conti';        % 'Micio Conti' or 'Chiusazza'

switch site
    case 'Micio Conti'
        % Held-out region: XY polygon
        polyXY = [ 510333.5  4157125;   % P1
                   510291.1  4157115;   % P2
                   510292.7  4157125;   % P3
                   510325.0  4157133 ]; % P4
        AOI.Zmin = 268.50;
        AOI.Zmax = 278.93;  % top of the cropped cloud (see 1_crop): the 914 points above 276.93 m
                            % belong to the held-out region (43,776 held-out / 21,966 constraint points)
        Zeval    = 276.93;  % top of the evaluation domain (top of the ERT model)
    case 'Chiusazza'
        % Held-out region: strip X 514160-514300, Y 4097700-4097730 (rectangle)
        % (20,510 held-out / 30,195 constraint points)
        polyXY = [ 514160  4097700;
                   514300  4097700;
                   514300  4097730;
                   514160  4097730 ];
        AOI.Zmin = 79.826;
        AOI.Zmax = 115.7055;
        Zeval    = AOI.Zmax;
    otherwise
        error('Unknown site: %s', site);
end
inputFile = [site '_cropped_ERTarea.txt'];

%% ---------- 2) READ DATA ----------
data = readmatrix(inputFile);
xyz  = data(:,1:3);

%% ---------- 3) HELD-OUT SELECTION ----------
insideXY = inpolygon(xyz(:,1), xyz(:,2), polyXY(:,1), polyXY(:,2));
insideZ  = xyz(:,3) >= AOI.Zmin & xyz(:,3) <= AOI.Zmax;
isHeldout = insideXY & insideZ;

constraintXYZ = xyz(~isHeldout,:);
heldoutXYZ    = xyz(isHeldout,:);

%% ---------- 4) STATISTICS ----------
fprintf('Total points: %d\n', size(xyz,1));
fprintf('CONSTRAINT points: %d (%.1f%%)\n', size(constraintXYZ,1), 100*size(constraintXYZ,1)/size(xyz,1));
fprintf('HELDOUT points: %d (%.1f%%)\n', size(heldoutXYZ,1), 100*size(heldoutXYZ,1)/size(xyz,1));

if size(heldoutXYZ,1) < 200
    warning('Very few points in the held-out region.');
end

%% ---------- 5) SAVE ----------
writematrix(constraintXYZ, [site '_constraint_only.txt'], 'Delimiter', ' ');
writematrix(heldoutXYZ,    [site '_heldout.txt'],         'Delimiter', ' ');
% Held-out region used in the evaluation. Micio Conti: 3D bounding box of the held-out
% points (Z 271.42-277.35 m). Chiusazza: nominal strip (Z range of the ERT model).
HELDOUT.Xmin = min(heldoutXYZ(:,1)); HELDOUT.Xmax = max(heldoutXYZ(:,1));
HELDOUT.Ymin = min(heldoutXYZ(:,2)); HELDOUT.Ymax = max(heldoutXYZ(:,2));
HELDOUT.Zmin = min(heldoutXYZ(:,3)); HELDOUT.Zmax = max(heldoutXYZ(:,3));
if strcmp(site,'Chiusazza')
    HELDOUT.Xmin = min(polyXY(:,1)); HELDOUT.Xmax = max(polyXY(:,1));
    HELDOUT.Ymin = min(polyXY(:,2)); HELDOUT.Ymax = max(polyXY(:,2));
    HELDOUT.Zmin = AOI.Zmin;         HELDOUT.Zmax = Zeval;
end
save([site '_holdout_region.mat'], 'polyXY', 'AOI', 'HELDOUT');

fprintf('\nFiles written: %s_constraint_only.txt (input of 3_2), %s_heldout.txt (reference), %s_holdout_region.mat\n', site, site, site);
