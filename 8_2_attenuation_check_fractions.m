%% 8_2_attenuation_check_fractions.m
% Extension of 8_1_attenuation_check.m. In addition to the statistics of 8_1, for C, H
% and U (same elements) it reports the number and fraction of elements of the
% constrained region with rho >= 25/50/75/100% of the imposed value.
% Output: attenuation_check_<site>.csv, attenuation_fractions_<site>.csv

clear; clc;

site = 'Chiusazza';        % 'Micio Conti' or 'Chiusazza'

switch site
    case 'Micio Conti'
        files = struct('U','Model_Micio Conti_U.txt', ...
                       'C','Model_Micio Conti_C_rho75000.txt', ...
                       'H','Model_Micio Conti_H_rho75000.txt');
        cloudC = 'Micio Conti_cropped_ERTarea.txt';          % full geometry
        cloudH = 'Micio Conti_constraint_only.txt';          % constraint region only
        alphaRadius = 1.35; rhoImposed = 75000;
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50; AOI.Zmax = 276.93;
    case 'Chiusazza'
        files = struct('U','Model_Chiusazza_U.txt', ...
                       'C','Model_Chiusazza_C_rho20000.txt', ...
                       'H','Model_Chiusazza_H_rho20000.txt');
        cloudC = 'Chiusazza_cropped_ERTarea.txt';
        cloudH = 'Chiusazza_constraint_only.txt';
        alphaRadius = 2.50; rhoImposed = 20000;
        polyXY = [ 514191.91  4097666.70;
                   514306.99  4097666.70;
                   514306.99  4097785.75;
                   514191.91  4097785.75 ];
        AOI.Zmin = 79.82; AOI.Zmax = 115.71;
end

rows = {}; frac = {};
for k = 1:2
    if k == 1, tag = 'C'; cloud = readmatrix(cloudC); else, tag = 'H'; cloud = readmatrix(cloudH); end
    cloud = cloud(:,1:3); cloud = cloud(all(isfinite(cloud),2),:);
    shp = alphaShape(cloud(:,1), cloud(:,2), cloud(:,3), alphaRadius);

    dU = readmatrix(files.U); dM = readmatrix(files.(tag));
    P = dM(:,1:3);
    keep = inpolygon(P(:,1),P(:,2),polyXY(:,1),polyXY(:,2)) & P(:,3)>=AOI.Zmin & P(:,3)<=AOI.Zmax;
    inside = false(size(P,1),1);
    inside(keep) = inShape(shp, P(keep,1), P(keep,2), P(keep,3));

    if size(dU,1) ~= size(dM,1) || max(abs(dU(:,1:3)-dM(:,1:3)),[],'all') > 1e-6
        warning('U and %s have different elements: the comparison with U is not valid', tag);
    end
    rM = dM(inside,4); rU = dU(inside,4);
    q = prctile(rM,[10 50 90]); qU = prctile(rU,[10 50 90]);
    fprintf('%s | %s: %d elements inside the constrained region\n', site, tag, sum(inside));
    fprintf('   model %s: rho (10/50/90 pct) = %.0f / %.0f / %.0f ; median/imposed = %.2f ; >= 50%% of imposed: %.1f%%\n', ...
        tag, q(1), q(2), q(3), q(2)/rhoImposed, 100*mean(rM >= 0.5*rhoImposed));
    fprintf('   U, same elements:  rho (10/50/90 pct) = %.0f / %.0f / %.0f\n', qU(1), qU(2), qU(3));
    rows(end+1,:) = {tag, sum(inside), q(1), q(2), q(3), q(2)/rhoImposed, qU(2)}; %#ok<SAGROW>
    for f = [0.25 0.5 0.75 1.0]
        nM = nnz(rM >= f*rhoImposed); nU = nnz(rU >= f*rhoImposed);
        fprintf('   >= %3.0f%% imposed: %s %d/%d (%.1f%%) | U %d/%d (%.1f%%)\n', 100*f, tag, nM, numel(rM), 100*nM/numel(rM), nU, numel(rU), 100*nU/numel(rU));
        frac(end+1,:) = {tag, 100*f, numel(rM), nM, 100*nM/numel(rM), nU, 100*nU/numel(rU)}; %#ok<SAGROW>
    end
end
Tf = cell2table(frac, 'VariableNames', {'model','pct_of_imposed','nElements','nModel','pctModel','nU','pctU'});
writetable(Tf, sprintf('attenuation_fractions_%s.csv', strrep(site,' ','')));
disp(Tf);

T = cell2table(rows, 'VariableNames', {'model','nElements','p10','median','p90','median_over_imposed','median_U'});
writetable(T, sprintf('attenuation_check_%s.csv', strrep(site,' ','')));
disp(T);
