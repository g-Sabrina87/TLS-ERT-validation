%% 7_1_bootstrap_UH_spatial_block_buffer.m
% Paired spatial block bootstrap, U vs H (switch `site`: Micio Conti or Chiusazza).
% Evaluation set: elements inside the held-out rectangle, farther than bufferRadius
% from the constraint points, within maxDistTLS of the full TLS cloud.
% Metrics: rank AUC (primary), IoU at a fixed threshold.
% Each run appends one row to a cumulative CSV (one row per block size).
clear; clc; close all;
site = 'Micio Conti';        % 'Micio Conti' or 'Chiusazza'
nBootstrap = 2000; rng(1);
blockSize = [3.0 3.0 2.0];   % m; paper: horizontal 3, 5, 7, 10 m (and 15, 20 m at Chiusazza), 2 m vertical
maxBootstrapAttempts = 100;

switch site
    case 'Micio Conti'
        siteName = 'Micio Conti H_75000';
        modelFile_U = 'Model_Micio Conti_U.txt';
        modelFile_H = 'Model_Micio Conti_H_rho75000.txt';
        tlsFullFile = 'Micio Conti_cropped_ERTarea.txt';
        constraintFile = 'Micio Conti_constraint_only.txt';
        holdoutFile = 'Micio Conti_holdout_region.mat';
        alphaRadius = 1.35;
        classificationThreshold = 48107;
        bufferRadius = 1.5;      % m (one mesh cell)
        polyXY = [ 510340.5  4157110;
                   510311.1  4157106;
                   510290.7  4157124;
                   510325.0  4157133 ];
        AOI.Zmin = 251.50;
        AOI.Zmax = 276.93;
        thresholdRangeAUC = linspace(1000, 60000, 50);
    case 'Chiusazza'
        siteName = 'Chiusazza H_20000';
        modelFile_U = 'Model_Chiusazza_U.txt';
        modelFile_H = 'Model_Chiusazza_H_rho20000.txt';
        tlsFullFile = 'Chiusazza_cropped_ERTarea.txt';
        constraintFile = 'Chiusazza_constraint_only.txt';
        holdoutFile = 'Chiusazza_holdout_region.mat';
        alphaRadius = 2.50;
        classificationThreshold = 5513;
        bufferRadius = 1.0;      % m (one mesh cell)
        polyXY = [ 514191.9125  4097666.704;
                   514306.9915  4097666.704;
                   514306.9915  4097785.755;
                   514191.9125  4097785.755 ];
        AOI.Zmin = 79.82597;
        AOI.Zmax = 115.70548;
        thresholdRangeAUC = linspace(1000, 20000, 50);
    otherwise
        error('Unknown site: %s', site);
end
load(holdoutFile, 'HELDOUT');
maxDistTLS = 3;              % m: keep only elements within this distance from the full TLS (Inf = no limit)

if ~isfile(modelFile_U)||~isfile(modelFile_H)||~isfile(tlsFullFile), error('Input file not found.'); end

tlsData=readmatrix(tlsFullFile); tlsXYZ=tlsData(:,1:3);
shpFull=alphaShape(tlsXYZ(:,1),tlsXYZ(:,2),tlsXYZ(:,3),alphaRadius);
[XYZ_HO_U,rho_HO_U,truth_HO_U]=extract_heldout(modelFile_U,polyXY,AOI,HELDOUT,shpFull);
[XYZ_HO_H,rho_HO_H,truth_HO_H]=extract_heldout(modelFile_H,polyXY,AOI,HELDOUT,shpFull);

if numel(rho_HO_U)~=numel(rho_HO_H), error('U and H have a different number of held-out elements.'); end
if size(XYZ_HO_U,1)~=size(XYZ_HO_H,1), error('U and H have different size.'); end
maxCoordDifference=max(abs(XYZ_HO_U(:)-XYZ_HO_H(:)));
if maxCoordDifference>1e-8, error('U/H coordinates do not match. Max diff=%.6g',maxCoordDifference); end
if any(truth_HO_U~=truth_HO_H), error('U/H TLS truth does not match.'); end

% --- EXCLUSION of the elements close to the CONSTRAINT points ---
cp = readmatrix(constraintFile); cp = cp(:,1:3);
dC = minDistToSet(XYZ_HO_U, cp);
dT = minDistToSet(XYZ_HO_U, tlsXYZ(:,1:3));
keepB = (dC > bufferRadius) & (dT <= maxDistTLS);
fprintf('Buffer %.1f m, within %.1f m of the TLS: elements kept %d out of %d\n', bufferRadius, maxDistTLS, sum(keepB), numel(keepB));
XYZ_HO_U = XYZ_HO_U(keepB,:); rho_HO_U = rho_HO_U(keepB); truth_HO_U = truth_HO_U(keepB);
XYZ_HO_H = XYZ_HO_H(keepB,:); rho_HO_H = rho_HO_H(keepB); truth_HO_H = truth_HO_H(keepB);

XYZ_HO=XYZ_HO_U; truth_HO=truth_HO_U; n=size(XYZ_HO,1);
blockID=make_spatial_blocks(XYZ_HO,blockSize);
[~,~,blockIndex]=unique(blockID); nBlocks=max(blockIndex);
blockCounts=accumarray(blockIndex,1);
fprintf('%s: n=%d, blocks=%d, block=[%.2f %.2f %.2f] m\n',siteName,n,nBlocks,blockSize);

% Rank-based AUC (Mann-Whitney): no threshold grid
classFixed_U=rho_HO_U>=classificationThreshold; classFixed_H=rho_HO_H>=classificationThreshold;

% Element indices of each block, computed once
blockElementIndices = accumarray(blockIndex, (1:n)', [nBlocks,1], @(v){v}, {[]});

IoU_U_boot=nan(nBootstrap,1); IoU_H_boot=nan(nBootstrap,1); AUC_U_boot=nan(nBootstrap,1); AUC_H_boot=nan(nBootstrap,1);
nValid=0; attempts=0;
while nValid<nBootstrap
    attempts=attempts+1; if attempts>nBootstrap*maxBootstrapAttempts, error('Unable to generate valid replicates. Reduce blockSize.'); end
    sampledBlocks=randi(nBlocks,[nBlocks 1]);
    idx=cell2mat(blockElementIndices(sampledBlocks));
    truth_b=truth_HO(idx); if all(truth_b)||~any(truth_b), continue; end
    nValid=nValid+1;
    IoU_U_boot(nValid)=compute_iou(classFixed_U(idx),truth_b); IoU_H_boot(nValid)=compute_iou(classFixed_H(idx),truth_b);
    AUC_U_boot(nValid)=auc_rank(rho_HO_U(idx),truth_b); AUC_H_boot(nValid)=auc_rank(rho_HO_H(idx),truth_b);
end

diffIoU=IoU_H_boot-IoU_U_boot; diffAUC=AUC_H_boot-AUC_U_boot;
ciIoU_U=prctile(IoU_U_boot,[2.5 97.5]); ciIoU_H=prctile(IoU_H_boot,[2.5 97.5]);
ciAUC_U=prctile(AUC_U_boot,[2.5 97.5]); ciAUC_H=prctile(AUC_H_boot,[2.5 97.5]);
ciDiffIoU=prctile(diffIoU,[2.5 97.5]); ciDiffAUC=prctile(diffAUC,[2.5 97.5]);

fprintf('\nIoU U=%.3f [%.3f, %.3f]\n',median(IoU_U_boot),ciIoU_U(1),ciIoU_U(2));
fprintf('IoU H=%.3f [%.3f, %.3f]\n',median(IoU_H_boot),ciIoU_H(1),ciIoU_H(2));
fprintf('Delta IoU H-U=%+.3f [%.3f, %.3f]\n',median(diffIoU),ciDiffIoU(1),ciDiffIoU(2));
fprintf('AUC U=%.3f [%.3f, %.3f]\n',median(AUC_U_boot),ciAUC_U(1),ciAUC_U(2));
fprintf('AUC H=%.3f [%.3f, %.3f]\n',median(AUC_H_boot),ciAUC_H(1),ciAUC_H(2));
fprintf('Delta AUC H-U=%+.3f [%.3f, %.3f]\n',median(diffAUC),ciDiffAUC(1),ciDiffAUC(2));

results=struct('siteName',siteName,'nElements',n,'nBlocks',nBlocks,'blockSize',blockSize,'classificationThreshold',classificationThreshold,'thresholdRangeAUC',thresholdRangeAUC,'IoU_U_boot',IoU_U_boot,'IoU_H_boot',IoU_H_boot,'AUC_U_boot',AUC_U_boot,'AUC_H_boot',AUC_H_boot,'diffIoU',diffIoU,'diffAUC',diffAUC,'ciIoU_U',ciIoU_U,'ciIoU_H',ciIoU_H,'ciAUC_U',ciAUC_U,'ciAUC_H',ciAUC_H,'ciDiffIoU',ciDiffIoU,'ciDiffAUC',ciDiffAUC,'blockCounts',blockCounts);
save(sprintf('bootstrap_spatial_block_%s_%gx%gx%g.mat',regexprep(siteName,'\s+','_'),blockSize(1),blockSize(2),blockSize(3)),'results');

% --- Cumulative CSV: one row per run (block size) ---
sensitivityFile = sprintf('spatial_block_sensitivity_buffer_%s_%gm_near%g.csv', regexprep(siteName,'\s+','_'), bufferRadius, maxDistTLS);
newRow = table(string(siteName), blockSize(1), blockSize(2), blockSize(3), nBlocks, n, ...
    median(diffIoU), ciDiffIoU(1), ciDiffIoU(2), median(diffAUC), ciDiffAUC(1), ciDiffAUC(2), ...
    'VariableNames', {'site','blockX','blockY','blockZ','nBlocks','nElements', ...
                       'medianDeltaIoU','ciDeltaIoU_low','ciDeltaIoU_high', ...
                       'medianDeltaAUC','ciDeltaAUC_low','ciDeltaAUC_high'});
if isfile(sensitivityFile)
    writetable(newRow, sensitivityFile, 'WriteMode', 'append');
else
    writetable(newRow, sensitivityFile);
end
fprintf('Row added to %s\n', sensitivityFile);

figure('Color','w'); subplot(1,2,1); histogram(diffIoU); hold on; xline(0,'k--'); xline(ciDiffIoU(1),'r:'); xline(ciDiffIoU(2),'r:'); xlabel('Delta IoU (H-U)'); grid on;
subplot(1,2,2); histogram(diffAUC); hold on; xline(0,'k--'); xline(ciDiffAUC(1),'r:'); xline(ciDiffAUC(2),'r:'); xlabel('Delta AUC (H-U)'); grid on;

function dmin = minDistToSet(P, Q)
o = mean(Q, 1); P = P - o; Q = Q - o; q2 = sum(Q.^2, 2)';
dmin = zeros(size(P,1), 1); chunk = 200;
for i = 1:chunk:size(P,1)
    ii = i:min(i+chunk-1, size(P,1)); Pc = P(ii,:);
    D2 = sum(Pc.^2, 2) + q2 - 2*(Pc*Q');
    dmin(ii) = sqrt(max(min(D2, [], 2), 0));
end
end
function [XYZ_HO,rho_HO,truth_HO]=extract_heldout(modelFile,polyXY,AOI,HELDOUT,shpFull)
data=readmatrix(modelFile); XYZ=data(:,1:3); rho=data(:,4);
inTOT=inpolygon(XYZ(:,1),XYZ(:,2),polyXY(:,1),polyXY(:,2)) & XYZ(:,3)>=AOI.Zmin & XYZ(:,3)<=AOI.Zmax;
XYZ_TOT=XYZ(inTOT,:); rho_TOT=rho(inTOT);
inHO=XYZ_TOT(:,1)>=HELDOUT.Xmin & XYZ_TOT(:,1)<=HELDOUT.Xmax & XYZ_TOT(:,2)>=HELDOUT.Ymin & XYZ_TOT(:,2)<=HELDOUT.Ymax & XYZ_TOT(:,3)>=HELDOUT.Zmin & XYZ_TOT(:,3)<=HELDOUT.Zmax;
XYZ_HO=XYZ_TOT(inHO,:); rho_HO=rho_TOT(inHO); truth_HO=inShape(shpFull,XYZ_HO(:,1),XYZ_HO(:,2),XYZ_HO(:,3));
end
function blockID=make_spatial_blocks(XYZ,blockSize)
XYZ0=min(XYZ,[],1); blockCoord=floor((XYZ-XYZ0)./blockSize); [~,~,blockID]=unique(blockCoord,'rows');
end
function iou=compute_iou(pred,truth)
intersection=sum(pred&truth); unionSet=sum(pred|truth); iou=intersection/max(unionSet,1);
end
function auc=auc_rank(score,truth)
% Exact AUC = P(score of a positive > score of a negative), ties = 1/2
nPos=sum(truth); nNeg=sum(~truth);
if nPos==0||nNeg==0, auc=NaN; return; end
[~,~,g]=unique(score); cnt=accumarray(g,1); avg=cumsum(cnt)-(cnt-1)/2; ranks=avg(g);
auc=(sum(ranks(truth))-nPos*(nPos+1)/2)/(nPos*nNeg);
end
