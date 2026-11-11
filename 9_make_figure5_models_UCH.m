%% 9_make_figure_models_UCH.m
% Figure B: plan-view maps (X-Y) of the models, one figure per site with three panels:
%   U (log10 rho) | C - U | H - U  (differences of log10 rho, common diverging scale).
% No interpolation: the model elements inside the area and the cavity depth range are
% drawn as squares (one tile per X-Y column, optionally smoothed for display, see site(k).smooth). Overlays: full TLS cloud (black), constraint cloud (blue), held-out
% sector (red rectangle) and evaluation-area polygon (dashed).
% The depth range is given by the 5th-95th percentiles of the full TLS cloud.
% Output: Figure_models_UCH_plan_MicioConti.png, Figure_models_UCH_plan_Chiusazza.png

clear; clc; close all;

site(1).name       = 'Micio Conti';
site(1).files      = {'Model_Micio Conti_U.txt', 'Model_Micio Conti_C_rho75000.txt', 'Model_Micio Conti_H_rho75000.txt'};   % U, C, H
site(1).tls        = 'Micio Conti_cropped_ERTarea.txt';
site(1).constraint = 'Micio Conti_constraint_only.txt';
site(1).holdout    = 'Micio Conti_holdout_region.mat';
site(1).polyXY     = [510340.5 4157110; 510311.1 4157106; 510290.7 4157124; 510325.0 4157133];
site(1).clim       = [3.4 4.8];  % log10(ohm*m): 5th-95th percentile of U in the plotted elements is 3.5-4.7
site(1).cell       = 1.471;      % mesh cell size (m): tiles drawn at true size; [] = scatter of squares
site(1).smooth     = true;       % smooth display of the tiles (linear interpolation on a 0.2 m grid); false = plain tiles
site(1).xmin       = [];
site(1).dlim       = 0.15;       % scale of the difference (log10 ohm*m); C-U 5th-95th pct -0.12..+0.25, H-U -0.05..+0.05

site(2).name       = 'Chiusazza';
site(2).files      = {'Model_Chiusazza_U.txt', 'Model_Chiusazza_C_rho20000.txt', 'Model_Chiusazza_H_rho20000.txt'};   % U, C, H
site(2).tls        = 'Chiusazza_cropped_ERTarea.txt';
site(2).constraint = 'Chiusazza_constraint_only.txt';
site(2).holdout    = 'Chiusazza_holdout_region.mat';
site(2).polyXY     = [514191.91 4097696.67; 514302.49 4097696.67; 514302.49 4097785.75; 514191.91 4097785.75];
site(2).clim       = [2.5 4.5];
site(2).cell       = 1.0;        % approximate cell size (m); only used for plain tiles
site(2).smooth     = true;
site(2).xmin       = min(site(2).polyXY(:,1));   % plots start at the western edge of the ERT area
site(2).dlim       = 0.3;

sitesToDo = [1 2];             % 1 = Micio Conti, 2 = Chiusazza
maxPts    = 40000;             % maximum number of TLS points drawn
titles    = {'U (unconstrained)', 'C - U (full geometry)', 'H - U (constraint region only)'};
t = linspace(0,1,256)'; divmap = [interp1([0 .5 1],[0.05 1 0.75],t) interp1([0 .5 1],[0.2 1 0.05],t) interp1([0 .5 1],[0.75 1 0.05],t)];   % blue-white-red, more saturated ends

zPct = [5 95];                 % depth range = percentiles of the full TLS (the cavity)
for s = sitesToDo
    S = site(s);
    tls = readmatrix(S.tls);        tls = tls(all(isfinite(tls(:,1:3)),2), 1:3);
    cp  = readmatrix(S.constraint); cp  = cp(all(isfinite(cp(:,1:3)),2), 1:3);
    L = load(S.holdout, 'HELDOUT'); HO = L.HELDOUT;
    zlo = prctile(tls(:,3), zPct(1)); zhi = prctile(tls(:,3), zPct(2));
    fprintf('%s: cavity depth range %.1f - %.1f m\n', S.name, zlo, zhi);

    dU = readmatrix(S.files{1}); PU = dU(:,1:3); rU = dU(:,4);
    inA = inpolygon(PU(:,1),PU(:,2),S.polyXY(:,1),S.polyXY(:,2)) & PU(:,3)>=zlo & PU(:,3)<=zhi;

    figure('Color','w','Position',[40 60 1650 520]);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    ix = round(linspace(1, size(tls,1), min(maxPts, size(tls,1))));
    jx = round(linspace(1, size(cp,1),  min(maxPts, size(cp,1))));
    for m = 1:3
        nexttile; hold on;
        d = readmatrix(S.files{m});
        if size(d,1) ~= size(dU,1) || max(abs(d(:,1:3)-PU),[],'all') > 1e-6, error('%s: elements different from U', S.files{m}); end
        if m == 1, val = log10(d(:,4)); else, val = log10(d(:,4)) - log10(rU); end
        ii = find(inA);
        if m == 1, [~,o] = sort(val(ii)); else, [~,o] = sort(abs(val(ii))); end   % larger values on top
        ii = ii(o);
        if ~isempty(S.cell)
            % regular mesh with several layers in the depth range: one tile per X-Y column,
            % drawn at the true cell size (mean of log10 rho for U; value of largest |change| for C-U, H-U)
            [~,~,g] = unique(round(PU(ii,1:2),1), 'rows');
            xc = accumarray(g, PU(ii,1), [], @mean);  yc = accumarray(g, PU(ii,2), [], @mean);
            if m == 1, vc = accumarray(g, val(ii), [], @mean);
            else,      vc = accumarray(g, val(ii), [], @(v) v(find(abs(v)==max(abs(v)),1))); end
            h = S.cell/2;
            if S.smooth
                % display only: column values interpolated linearly on a fine grid and clipped to the area
                % (the model itself is not interpolated; the tiles above are the actual values)
                F = scatteredInterpolant(xc, yc, vc, 'linear', 'none');
                [Xg, Yg] = meshgrid(min(xc):0.2:max(xc), min(yc):0.2:max(yc));
                Vg = F(Xg, Yg);
                Vg(~inpolygon(Xg, Yg, S.polyXY(:,1), S.polyXY(:,2))) = NaN;
                pcolor(Xg, Yg, Vg); shading interp;
            else
                patch((xc + h*[-1 1 1 -1])', (yc + h*[-1 -1 1 1])', vc', 'EdgeColor', 'none');
            end
        else
            scatter(PU(ii,1), PU(ii,2), 14, val(ii), 's', 'filled');
        end
        plot(tls(ix,1), tls(ix,2), '.', 'Color', [0 0 0], 'MarkerSize', 1);
        plot(cp(jx,1),  cp(jx,2),  '.', 'Color', [0.1 0.3 0.9], 'MarkerSize', 1);
        plot([HO.Xmin HO.Xmax HO.Xmax HO.Xmin HO.Xmin], [HO.Ymin HO.Ymin HO.Ymax HO.Ymax HO.Ymin], '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 2);
        plot([S.polyXY(:,1); S.polyXY(1,1)], [S.polyXY(:,2); S.polyXY(1,2)], 'k--', 'LineWidth', 0.8);
        axis equal tight; box on; xlabel('X (m)'); ylabel('Y (m)');
        if ~isempty(S.xmin), xl = xlim; xlim([S.xmin xl(2)]); end
        ax = gca; ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0; xtickformat('%.0f'); ytickformat('%.0f');   % full UTM coordinates, no x10^n
        ax.XTickLabelRotation = 30;
        if m == 1, caxis(S.clim); colormap(gca, parula); cb1 = colorbar('southoutside'); cb1.Label.String = 'log_{10} \rho (\Omega\cdotm)';
        else, caxis([-S.dlim S.dlim]); colormap(gca, divmap); cb2 = colorbar('southoutside'); cb2.Label.String = sprintf('\\Delta log_{10} \\rho vs U (saturated at \\pm%.2f)', S.dlim); end
        title(titles{m}, 'FontWeight', 'normal');
    end
    title(tl, sprintf('%s - plan view, elements at Z = %.0f-%.0f m; red box = held-out sector, blue = constraint region, black = TLS', S.name, zlo, zhi), 'FontSize', 10);
    outName = sprintf('Figure_models_UCH_plan_%s.png', strrep(S.name, ' ', ''));
    exportgraphics(gcf, outName, 'Resolution', 300);
    fprintf('%s written\n', outName);
end
