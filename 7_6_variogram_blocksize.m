%% 7_6_variogram_blocksize.m
% Spatial autocorrelation scale of log-resistivity of model U on the evaluation set,
% to compare with the bootstrap block sizes (3, 5, 7, 10 m). Experimental variogram
% with a spherical model (nugget, sill, range) on random subsamples; two versions:
% isotropic 3D and horizontal (|dz| <= zTol).
% Output: variogram_<site>.csv and variogram_<site>.png

clear; clc;

site = 'Chiusazza';        % 'Micio Conti' or 'Chiusazza'
nSub = 3000;               % points per subsample
nRep = 5;                  % independent subsamples
zTol = 1.5;                % m, vertical tolerance for the horizontal version
lagEdges = 0:1:30;         % m

switch site
    case 'Chiusazza'
        modelFile_U = 'Model_Chiusazza_U.txt';
        holdoutFile = 'Chiusazza_holdout_region.mat';
        constraintFile = 'Chiusazza_constraint_only.txt';
        tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';
        polyXY = [ 514191.91  4097666.70;
                   514306.99  4097666.70;
                   514306.99  4097785.75;
                   514191.91  4097785.75 ];
        AOI.Zmin = 79.82;
        AOI.Zmax = 115.71; bufferRadius = 1.0;
     case 'Micio Conti'
        modelFile_U = 'Model_Micio Conti_U.txt';
        holdoutFile = 'Micio Conti_holdout_region.mat';
        constraintFile = 'Micio Conti_constraint_only.txt';
        tlsFullFile = 'Micio Conti_cropped_ERTarea.txt';
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93; bufferRadius = 1.5;
end
maxDistTLS = 3;
rng(1);

%% ---------- EVALUATION SET (as in the paper) ----------
load(holdoutFile,'HELDOUT');
cp  = readmatrix(constraintFile); cp  = cp(all(isfinite(cp(:,1:3)),2),1:3);
tls = readmatrix(tlsFullFile);    tls = tls(all(isfinite(tls(:,1:3)),2),1:3);
d = readmatrix(modelFile_U); XYZ = d(:,1:3); lr = log10(d(:,4));
inA = inpolygon(XYZ(:,1),XYZ(:,2),polyXY(:,1),polyXY(:,2)) & XYZ(:,3)>=AOI.Zmin & XYZ(:,3)<=AOI.Zmax;
inHO = XYZ(:,1)>=HELDOUT.Xmin & XYZ(:,1)<=HELDOUT.Xmax & XYZ(:,2)>=HELDOUT.Ymin & XYZ(:,2)<=HELDOUT.Ymax & ...
       XYZ(:,3)>=HELDOUT.Zmin & XYZ(:,3)<=HELDOUT.Zmax;
rect = inA & inHO; P = XYZ(rect,:); v = lr(rect);
[~,dC] = knnsearch(cp,P,'K',1); [~,dT] = knnsearch(tls,P,'K',1);
ev = (dC > bufferRadius) & (dT <= maxDistTLS);
P = P(ev,:); v = v(ev);
fprintf('%s: evaluation set = %d elements (expected: 355 Micio Conti, 11758 Chiusazza)\n', site, size(P,1));

%% ---------- VARIOGRAM ----------
nL = numel(lagEdges)-1; lagMid = (lagEdges(1:end-1)+lagEdges(2:end))/2;
G3 = nan(nRep,nL); Gh = nan(nRep,nL); N3 = zeros(nRep,nL); Nh = zeros(nRep,nL);
for r = 1:nRep
    ii = randperm(size(P,1), min(nSub,size(P,1)));
    Ps = P(ii,:); vs = v(ii);
    h3 = pdist(Ps); sq = pdist(vs,'squaredeuclidean')/2;
    hxy = pdist(Ps(:,1:2)); dz = pdist(Ps(:,3));
    for b = 1:nL
        s = h3 >= lagEdges(b) & h3 < lagEdges(b+1);
        N3(r,b) = nnz(s); if N3(r,b) > 0, G3(r,b) = mean(sq(s)); end
        s = hxy >= lagEdges(b) & hxy < lagEdges(b+1) & dz <= zTol;
        Nh(r,b) = nnz(s); if Nh(r,b) > 0, Gh(r,b) = mean(sq(s)); end
    end
end

%% ---------- SPHERICAL FIT ----------
sph = @(p,h) p(1) + (p(2)-p(1)) .* ( (h<p(3)).*(1.5*h./p(3) - 0.5*(h./p(3)).^3) + (h>=p(3)) );
fit1 = @(g,n,tag) fitOne(g,n,lagMid,sph,tag,var(v));
res3 = nan(nRep,3); resh = nan(nRep,3);
for r = 1:nRep
    res3(r,:) = fit1(G3(r,:),N3(r,:),'3D');
    resh(r,:) = fit1(Gh(r,:),Nh(r,:),'horiz.');
end
fprintf('\nRange of the spherical model (m), %d subsamples:\n', nRep);
fprintf('  isotropic 3D : median %.1f  (min %.1f, max %.1f)\n', median(res3(:,3)), min(res3(:,3)), max(res3(:,3)));
fprintf('  horizontal   : median %.1f  (min %.1f, max %.1f)\n', median(resh(:,3)), min(resh(:,3)), max(resh(:,3)));
fprintf('Sample variance of log10(rho): %.4f ; fitted sill (3D): %.4f ; nugget: %.4f\n', var(v), median(res3(:,2)), median(res3(:,1)));

T = table(lagMid', mean(G3,1,'omitnan')', sum(N3,1)', mean(Gh,1,'omitnan')', sum(Nh,1)', ...
    'VariableNames',{'lag_m','gamma_3D','pairs_3D','gamma_horiz','pairs_horiz'});
writetable(T, sprintf('variogram_%s.csv', strrep(site,' ','')));
figure; hold on; plot(lagMid, mean(G3,1,'omitnan'), 'o-'); plot(lagMid, mean(Gh,1,'omitnan'), 's-');
xline(median(res3(:,3)),'--'); yline(var(v),':');
xlabel('Lag (m)'); ylabel('Semivariance of log_{10}\rho'); legend('3D','horizontal','3D range','variance'); title(site); grid on;
saveas(gcf, sprintf('variogram_%s.png', strrep(site,' ','')));

%% ---------- LOCAL FUNCTION ----------
function p = fitOne(g,n,h,sph,~,v0)
    ok = isfinite(g) & n > 50;
    if nnz(ok) < 5, p = [NaN NaN NaN]; return; end
    w = sqrt(n(ok));
    cost = @(q) sum( w.^2 .* (sph([abs(q(1)) abs(q(2)) abs(q(3))+0.5], h(ok)) - g(ok)).^2 );
    q0 = [0.2*v0, v0, 8];
    q = fminsearch(cost, q0, optimset('Display','off','MaxIter',2000));
    p = [abs(q(1)) abs(q(2)) abs(q(3))+0.5];
end
