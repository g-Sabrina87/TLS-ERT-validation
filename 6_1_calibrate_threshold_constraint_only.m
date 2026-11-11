%% 6_1_calibrate_threshold_constraint_only.m
% Calibrates the classification threshold used in the U/C/H comparison with no
% leakage from the held-out zone:
%   - the reference volume comes ONLY from the CONSTRAINT geometry (not from the
%     full TLS, which would include the held-out zone);
%   - the resistivity field on which the threshold is searched comes from the
%     UNCONSTRAINED model (U), not from C or H.
% The threshold is the resistivity whose isosurface encloses the same volume as the
% alpha shape of the constraint cloud.

clear; clc; close all;

%% ---------- 1) PARAMETERS ----------
siteName = 'Micio Conti';        % 'Micio Conti' or 'Chiusazza'

switch siteName
    case 'Micio Conti'
        modelFile_U = 'Model_Micio Conti_U.txt';                 % unconstrained model
        constraintTlsFile = 'Micio Conti_constraint_only.txt';   % constraint geometry only (no held-out)
        alphaRadius = 1.35;
        % ERT area (polygon + Z range)
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93;
    case 'Chiusazza'
        modelFile_U = 'Model_Chiusazza_U.txt';
        constraintTlsFile = 'Chiusazza_constraint_only.txt';
        alphaRadius = 2.50;
        polyXY = [ 514191.9125  4097666.704;
                   514306.9915  4097666.704;
                   514306.9915  4097785.755;
                   514191.9125  4097785.755 ];
        AOI.Zmin = 79.82597;
        AOI.Zmax = 115.70548;
    otherwise
        error('Unknown site: %s', siteName);
end

gridStep = 1.0;

%% ---------- 2) REFERENCE VOLUME (constraint geometry only) ----------
constraintData = readmatrix(constraintTlsFile);
constraintXYZ  = constraintData(:,1:3);
shpConstraint  = alphaShape(constraintXYZ(:,1), constraintXYZ(:,2), constraintXYZ(:,3), alphaRadius);
volConstraint  = volume(shpConstraint);
fprintf('%s - volume of the CONSTRAINT geometry (reference): %.1f m^3\n', siteName, volConstraint);

%% ---------- 3) RESISTIVITY FIELD FROM THE UNCONSTRAINED MODEL (U) ----------
dataU = readmatrix(modelFile_U);
XYZ_U = dataU(:,1:3);
rho_U = dataU(:,4);

inAOI = inpolygon(XYZ_U(:,1), XYZ_U(:,2), polyXY(:,1), polyXY(:,2)) & ...
        XYZ_U(:,3) >= AOI.Zmin & XYZ_U(:,3) <= AOI.Zmax;
XYZ_U = XYZ_U(inAOI,:); rho_U = rho_U(inAOI);

xv = min(XYZ_U(:,1)):gridStep:max(XYZ_U(:,1));
yv = min(XYZ_U(:,2)):gridStep:max(XYZ_U(:,2));
zv = min(XYZ_U(:,3)):gridStep:max(XYZ_U(:,3));
[Xg, Yg, Zg] = meshgrid(xv, yv, zv);

Frho = scatteredInterpolant(XYZ_U(:,1), XYZ_U(:,2), XYZ_U(:,3), rho_U, 'linear', 'none');
Rhog = Frho(Xg, Yg, Zg);

%% ---------- 4) THRESHOLD SEARCH (volume matching with the constraint only) ----------
cellVolume = gridStep^3;
lowBound  = prctile(rho_U, 1);
highBound = prctile(rho_U, 99.9);
rhoRange  = linspace(lowBound, highBound, 300);

testVolume = nan(size(rhoRange));
for i = 1:numel(rhoRange)
    testVolume(i) = sum(Rhog(:) >= rhoRange(i), 'omitnan') * cellVolume;
end

[diffMin, idxMin] = min(abs(testVolume - volConstraint));
calibratedThreshold = rhoRange(idxMin);

fprintf('\n=== CALIBRATED THRESHOLD - %s ===\n', siteName);
fprintf('Threshold = %.1f ohm-m (isosurface volume = %.1f m^3, constraint target = %.1f m^3, diff=%.1f%%)\n', ...
    calibratedThreshold, testVolume(idxMin), volConstraint, 100*diffMin/volConstraint);

%% ---------- 5) PLOT ----------
figure;
plot(rhoRange, testVolume, 'o-'); hold on;
yline(volConstraint, 'r--', 'Constraint volume (target)');
plot(calibratedThreshold, testVolume(idxMin), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
xlabel('Tested threshold (ohm\cdotm)'); ylabel('Enclosed volume (m^3)');
title(sprintf('%s: threshold calibration (threshold = %.0f)', siteName, calibratedThreshold));
grid on;
