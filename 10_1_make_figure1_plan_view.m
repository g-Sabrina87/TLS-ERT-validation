%% 10_1_make_figureA_plan_view.m
% Figure A of the paper: plan view of each site with the ERT electrodes, the TLS point
% cloud, the constraint region (blue), the held-out sector (red box), the investigated
% area used for the model (dashed) and the evaluation elements of the main analysis.
% One figure per site; set sitesToDo = [1 2] to produce both.
% Evaluation elements (Section 3.5 of the paper): element centroids inside the held-out
% region, the investigated area and the depth range, farther than 'bufferRadius' from every
% constraint point and within 'maxDistTLS' of the full TLS cloud.
% Requires the Statistics and Machine Learning Toolbox (knnsearch).
% Output: Figure_plan_view_MicioConti.png and Figure_plan_view_Chiusazza.png

clear; clc; close all;

site(1).name       = 'Micio Conti';
site(1).tls        = 'Micio Conti_cropped_ERTarea.txt';                 % full cropped TLS cloud
site(1).constraint = 'Micio Conti_constraint_only.txt';                 % constraint-region cloud
site(1).model      = 'Model_Micio Conti_U.txt';                      % model U (element centroids: X Y Z rho)
site(1).electrodes = 'electrodes_MicioConti.csv';                  % cable,id,X,Y,Z
site(1).holdout    = 'Micio Conti_holdout_region.mat';                      % struct HELDOUT (Xmin Xmax Ymin Ymax Zmin Zmax)
site(1).polyXY     = [510340.5 4157110; 510311.1 4157106; 510290.7 4157124; 510325.0 4157133];
site(1).Z          = [251.50 276.93];
site(1).buffer     = 1.5;
site(1).tlsOnTop   = false;
site(1).legLoc     = 'northeast';   % legend position inside the plot (change if it covers data)

site(2).name       = 'Chiusazza';
site(2).tls        = 'Chiusazza_cropped_ERTarea.txt';
site(2).constraint = 'Chiusazza_constraint_only.txt';
site(2).model      = 'Model_Chiusazza_U.txt';
site(2).electrodes = 'electrodes_Chiusazza.csv';
site(2).holdout    = 'Chiusazza_holdout_region.mat';
site(2).polyXY     = [514191.91 4097696.67; 514302.49 4097696.67; 514302.49 4097785.75; 514191.91 4097785.75];
site(2).Z          = [79.826 115.7055];
site(2).buffer     = 1.0;
site(2).tlsOnTop   = true;
site(2).legLoc     = 'northeast';

sitesToDo   = [1 2];
maxDistTLS  = 3;        % m
maxPlotPts  = 30000;    % maximum number of TLS points drawn per cloud

for s = sitesToDo
    S = site(s);
    tls = readmatrix(S.tls);        tls = tls(all(isfinite(tls(:,1:3)),2), 1:3);
    cp  = readmatrix(S.constraint); cp  = cp(all(isfinite(cp(:,1:3)),2), 1:3);
    El  = readmatrix(S.electrodes); El  = El(:,3:5);
    L   = load(S.holdout, 'HELDOUT'); HO = L.HELDOUT;
    dU  = readmatrix(S.model); P = dU(:,1:3);

    % evaluation elements
    inA  = inpolygon(P(:,1), P(:,2), S.polyXY(:,1), S.polyXY(:,2)) & P(:,3) >= S.Z(1) & P(:,3) <= S.Z(2);
    inHO = P(:,1) >= HO.Xmin & P(:,1) <= HO.Xmax & P(:,2) >= HO.Ymin & P(:,2) <= HO.Ymax & ...
           P(:,3) >= HO.Zmin & P(:,3) <= HO.Zmax;
    cand = find(inA & inHO);
    [~, dC] = knnsearch(cp,  P(cand,:), 'K', 1);
    [~, dT] = knnsearch(tls, P(cand,:), 'K', 1);
    ev = cand(dC > S.buffer & dT <= maxDistTLS);
    fprintf('%s: %d evaluation elements (check against Section 3.5 of the paper)\n', S.name, numel(ev));

    ix = round(linspace(1, size(tls,1), min(maxPlotPts, size(tls,1))));
    jx = round(linspace(1, size(cp,1),  min(maxPlotPts, size(cp,1))));

    % held-out rectangle clipped to the investigated area (display only; no TLS point or
    % element lies outside the investigated area, so the evaluation is not affected)
    xr = [max(HO.Xmin, min(S.polyXY(:,1))) min(HO.Xmax, max(S.polyXY(:,1)))];
    yr = [max(HO.Ymin, min(S.polyXY(:,2))) min(HO.Ymax, max(S.polyXY(:,2)))];
    isHOtls = tls(:,1) >= HO.Xmin & tls(:,1) <= HO.Xmax & tls(:,2) >= HO.Ymin & tls(:,2) <= HO.Ymax & ...
              tls(:,3) >= HO.Zmin & tls(:,3) <= HO.Zmax;
    th = tls(isHOtls,:);  thx = round(linspace(1, size(th,1), min(maxPlotPts, size(th,1))));

    figure('Color', 'w', 'Position', [60 60 900 800]); hold on; box on;
    % evaluation elements are 3-D: in plan view the elements of all depths are superimposed,
    % hence semi-transparent markers (darker = more elements above each other)
    plot(tls(ix,1), tls(ix,2), '.', 'Color', [0.65 0.65 0.65], 'MarkerSize', 3);
    plot(cp(jx,1),  cp(jx,2),  '.', 'Color', [0.10 0.30 0.90], 'MarkerSize', 3);
    plot(th(thx,1), th(thx,2), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 3);
    scatter(P(ev,1), P(ev,2), 22, [0.95 0.55 0.0], 's', 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', [0.75 0.35 0.0], 'MarkerEdgeAlpha', 0.5, 'LineWidth', 0.3);   % drawn above the TLS clouds
    if S.tlsOnTop, plot(th(thx,1), th(thx,2), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 2); end   % keep the held-out TLS points visible above the elements
    plot(El(:,1),   El(:,2),   'k^', 'MarkerSize', 4, 'MarkerFaceColor', 'w');
    plot([xr(1) xr(2) xr(2) xr(1) xr(1)], [yr(1) yr(1) yr(2) yr(2) yr(1)], '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2);
    plot([S.polyXY(:,1); S.polyXY(1,1)], [S.polyXY(:,2); S.polyXY(1,2)], 'k--', 'LineWidth', 1);
    axis equal tight; xlabel('Easting, UTM 33N (m)'); ylabel('Northing, UTM 33N (m)');
    ax = gca; ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0; ax.XAxis.TickLabelFormat = '%.0f'; ax.YAxis.TickLabelFormat = '%.0f';
    xl = xlim; ax.XTick = ceil(xl(1)/20)*20 : 20 : xl(2);   % tick every 20 m (avoids overlapping labels)
    title(S.name, 'FontWeight', 'normal');
    % legend with dummy handles (larger markers than the data)
    d(1) = plot(NaN, NaN, '.', 'Color', [0.10 0.30 0.90], 'MarkerSize', 18);
    d(2) = plot(NaN, NaN, '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 18);
    d(3) = plot(NaN, NaN, 's', 'MarkerFaceColor', [0.95 0.55 0.0], 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    d(4) = plot(NaN, NaN, 'k^', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    d(5) = plot(NaN, NaN, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2);
    d(6) = plot(NaN, NaN, 'k--', 'LineWidth', 1);
    legend(d, {'Constraint region', 'TLS points of the held-out sector', ...
               'Evaluation elements (all depths)', 'ERT electrodes', 'Held-out sector', 'Evaluation area'}, 'Location', S.legLoc, 'FontSize', 8, 'Color', 'w', 'EdgeColor', [0.5 0.5 0.5]);
    exportgraphics(gcf, sprintf('Figure_plan_view_%s.png', strrep(S.name, ' ', '')), 'Resolution', 300);
end
