%% 3_3_1_chiusazza_extra_holdouts.m
% Additional held-out sectors for Chiusazza: strips of fixed width over the whole
% width of the domain, along Y.
%
% Two modes (variable 'modo'):
%   'lista'  -> for each candidate strip (sliding window along Y) computes, without
%               inversions, the TLS points in the strip, the constraint share and the
%               number of evaluation elements (exclusion 1.0 m from the constraint
%               points, within 3 m of the TLS). Row 1 is the current strip (check:
%               20,510 held-out points, 30,195 constraint points, 11,758 elements,
%               2,200 in the cavity).
%   'scrivi' -> for each strip in 'scelte' creates a folder with
%                 grotta_constraint_only_Chiusazza.txt  (constraint cloud)
%                 grotta_heldout_Chiusazza.txt          (held-out cloud, reference)
%                 holdout_region.mat and holdout_region_chiusazza.mat (same content)
%               holdout_region.mat is read by 3_2_4_2_cavity_to_rect_anomalies.m, to be
%               run inside the folder of the strip.
%
% Selection rule (set before any inversion): strips not overlapping the current strip
% nor each other, with at least 1,000 cavity elements in the evaluation set. 30 m
% strips were not applicable (only 55.75 m remain north of the current strip), so two
% adjacent strips of 27.875 m covering the domain north of the current strip were used.

clear; clc;

%% ---------- SETTINGS ----------
modo = 'scrivi';                                % 'lista' or 'scrivi'
tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';     % full TLS cloud (cropped to the ERT area)
modelFile_U = 'Model_Chiusazza_U.txt';           % model file (only the element coordinates are used)
alphaRadius = 2.50;
polyXY = [ 514191.91  4097696.67;
           514302.49  4097696.67;
           514302.49  4097785.75;
           514191.91  4097785.75 ];
AOI.Zmin = 79.8260; AOI.Zmax = 115.7055;        % Z range of the evaluation elements
HOZ = [79.8260 115.7055];                       % Z range of the held-out TLS points
HOX = [514160 514310];                          % X range of the new strips
holdoutFile = 'Chiusazza_holdout_region.mat';   % current strip (only for the check row)
bufferRadius = 1.0;                             % m, exclusion around the constraint points
maxDistTLS   = 3;                               % m, limit from the full TLS cloud
stripWidth   = 30;                              % m
stepY        = 5;                               % m, window step
YrefA        = [4097700 4097730];               % current strip (check)
minCavElems  = 1000;                            % selection rule

% Strips to write ([Ymin Ymax], one row per strip)
scelte = [ 4097730      4097757.875;
           4097757.875  4097785.75 ];

%% ---------- DATA ----------
tls = readmatrix(tlsFullFile); tls = tls(all(isfinite(tls(:,1:3)),2),1:3);
d = readmatrix(modelFile_U); XYZ = d(:,1:3);
inA = inpolygon(XYZ(:,1),XYZ(:,2),polyXY(:,1),polyXY(:,2)) & XYZ(:,3)>=AOI.Zmin & XYZ(:,3)<=AOI.Zmax;
XYZ = XYZ(inA,:);
Ydom = [min(polyXY(:,2)) max(polyXY(:,2))];
shp = alphaShape(tls(:,1),tls(:,2),tls(:,3),alphaRadius);
inCavAll = inShape(shp, XYZ(:,1), XYZ(:,2), XYZ(:,3));
fprintf('TLS points: %d ; model elements inside the AOI: %d\n', size(tls,1), size(XYZ,1));

%% ---------- STATISTICS OF A STRIP ----------
mkBox = @(Y0,Y1) struct('Xmin',HOX(1),'Xmax',HOX(2),'Ymin',Y0,'Ymax',Y1,'Zmin',HOZ(1),'Zmax',HOZ(2));
stat = @(H) stripStats(H,tls,XYZ,inCavAll,AOI,bufferRadius,maxDistTLS);
L = load(holdoutFile,'HELDOUT');   % current strip: X 514160-514300, Y 4097700-4097730, Z 79.826-115.7055

switch modo
    case 'lista'
        starts = Ydom(1):stepY:(Ydom(2)-stripWidth);
        rows = [YrefA, stat(L.HELDOUT)];
        for y0 = starts
            rows(end+1,:) = [y0, y0+stripWidth, stat(mkBox(y0,y0+stripWidth))]; %#ok<SAGROW>
        end
        T = array2table(rows,'VariableNames',{'Ymin','Ymax','ptsHeldout','ptsConstraint','pctConstraint','ptsOutsideZ','nEval','nCavity','pctCavity','overlapsA'});
        T.overlapsA = double(T.Ymin < YrefA(2) & T.Ymax > YrefA(1));
        T.overlapsA(1) = 0;
        writetable(T,'candidate_strips_Chiusazza.csv');
        disp(T);
        fprintf('\nRow 1 = current strip (check: 20510 / 30195 / 11758 / 2200).\n');
        ok = T(2:end,:); ok = ok(ok.overlapsA==0 & ok.nCavity>=minCavElems,:);
        fprintf('Candidates not overlapping the current strip with >= %d cavity elements: %d\n', minCavElems, height(ok));
        disp(sortrows(ok,'nCavity','descend'));
        fprintf('\nStrips in ''scelte'':\n');
        R2 = [];
        for i = 1:size(scelte,1)
            R2(end+1,:) = [scelte(i,:), stat(mkBox(scelte(i,1),scelte(i,2)))]; %#ok<SAGROW>
        end
        T2 = array2table(R2,'VariableNames',{'Ymin','Ymax','ptsHeldout','ptsConstraint','pctConstraint','ptsOutsideZ','nEval','nCavity','pctCavity','x'});
        T2.x = [];
        disp(T2);

    case 'scrivi'
        for i = 1:size(scelte,1)
            Y0 = scelte(i,1); Y1 = scelte(i,2);
            tag = sprintf('strip_Y%d_%d', round(Y0), round(Y1));
            if ~exist(tag,'dir'), mkdir(tag); end
            HELDOUT = mkBox(Y0,Y1);
            isHO = tls(:,1)>=HELDOUT.Xmin & tls(:,1)<=HELDOUT.Xmax & tls(:,2)>=Y0 & tls(:,2)<=Y1 & tls(:,3)>=HELDOUT.Zmin & tls(:,3)<=HELDOUT.Zmax;
            writematrix(tls(~isHO,:),fullfile(tag,'grotta_constraint_only_Chiusazza.txt'),'Delimiter',' ');
            writematrix(tls(isHO,:), fullfile(tag,'grotta_heldout_Chiusazza.txt'),'Delimiter',' ');
            % holdout_region.mat: polyXY and AOI describe the HELDOUT region (rectangle + Z)
            polyXY_ERT = polyXY; polyXY = [HELDOUT.Xmin HELDOUT.Ymin; HELDOUT.Xmax HELDOUT.Ymin; HELDOUT.Xmax HELDOUT.Ymax; HELDOUT.Xmin HELDOUT.Ymax];
            AOI_ERT = AOI; AOI.Zmin = HELDOUT.Zmin; AOI.Zmax = HELDOUT.Zmax;
            save(fullfile(tag,'holdout_region.mat'),'HELDOUT','polyXY','AOI');
            save(fullfile(tag,'holdout_region_chiusazza.mat'),'HELDOUT','polyXY','AOI');
            polyXY = polyXY_ERT; AOI = AOI_ERT;
            s = stat(HELDOUT);
            fprintf('%s: held-out %d points, constraint %d (%.1f%%), evaluation elements %d (%d in the cavity)\n', ...
                tag, s(1), s(2), s(3), s(5), s(6));
            if s(4) > 0, warning('%s: %d TLS points of the strip are outside the Z range: they remain in the constraint.', tag, s(4)); end
        end
        fprintf('\nNext, in each folder: 3_2_4_2_cavity_to_rect_anomalies.m (alphaRadius 2.5, voxelSize 1.0, rho 20000)\n');
        fprintf('-> ERTLab boxes -> H inversion -> evaluation.\n');
end

%% ---------- LOCAL FUNCTION ----------
function out = stripStats(H,tls,XYZ,inCavAll,AOI,bufferRadius,maxDistTLS)
    inStripXY = tls(:,1)>=H.Xmin & tls(:,1)<=H.Xmax & tls(:,2)>=H.Ymin & tls(:,2)<=H.Ymax;
    isHO = inStripXY & tls(:,3)>=H.Zmin & tls(:,3)<=H.Zmax;
    ptsOutsideZ = nnz(inStripXY & ~isHO);          % points in the strip (XY) outside the Z range: they would remain in the constraint
    cp = tls(~isHO,:);
    rect = XYZ(:,1)>=H.Xmin & XYZ(:,1)<=H.Xmax & XYZ(:,2)>=H.Ymin & XYZ(:,2)<=H.Ymax;   % elements already limited to the AOI (polygon and Z)
    P = XYZ(rect,:);
    if isempty(cp) || isempty(P)
        out = [nnz(isHO) size(cp,1) 100*size(cp,1)/size(tls,1) ptsOutsideZ 0 0 0 0]; return;
    end
    [~,dC] = knnsearch(cp, P, 'K', 1);
    [~,dT] = knnsearch(tls, P, 'K', 1);
    ev = (dC > bufferRadius) & (dT <= maxDistTLS);
    cav = inCavAll(rect); cav = cav(ev);
    out = [nnz(isHO), size(cp,1), 100*size(cp,1)/size(tls,1), ptsOutsideZ, nnz(ev), nnz(cav), 100*nnz(cav)/max(nnz(ev),1), 0];
end
