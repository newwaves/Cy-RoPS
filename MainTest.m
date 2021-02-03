%-------------------------------------------------------------------------------
%   Main function for localization
%   1-Delete ground --> 2-Feature --> 3-Delete Obstacles --> localization
%   modify by Z.X.Tao E_mail: newwaves@stu.xjtu.edu.cn-- 2018-02-01 
%-------------------------------------------------------------------------------
clc; close all;clear all;
FixPt = [300,300,-20;...
        300,300,20;...
        300,-300,-20;...
        300,-300,20;...
        -300,300,-20;...
        -300,300,20;...
        -300,-300,-20;...
        -300,-300,20];
MapName = '.\Map\GarageMap.pcd';  %   GarageWithNormal_0.2.ply
ptCloud_Map = pcread( MapName );             
figure;hold on;grid on;
pcshow(ptCloud_Map);
xlabel('x');ylabel('y');zlabel('z');
title('Map');
ptCloudB_Reduced = pcdownsample(ptCloud_Map, 'random', 0.5); % for display
KdMap = createns(ptCloud_Map.Location);
%% --------------------------------------------------------------------
isShow = 0;
load '.\Cy-RoPS\CylinderFv.mat';
CylinderFv = mean(cc);
load '.\Cy-RoPS\WallFv.mat';
WallFv = mean(cc);
load '.\Cy-RoPS\VehicleFv.mat';
VehicleFv = mean(cc);
%%
HDLFolder = '.\Data';
Len = length(dir(HDLFolder)) * 10;
% ----------------------------
AssumTime = []; % zeros( FrameLength,1 );
vRmse = [];
TRatio = [];
vPtNum = [];
% --------------------------------------------------------------------
TraQuat = [-25.1443,38.3366,-4.57785,-0.00326366,0.00395295,-0.751049,0.660226];% Initial pose
TTf = CQuat2TraTTf(TraQuat); 
R = eul2rotm(rotm2eul(TTf(1:3, 1:3)));
T = TTf(1:3, end);
TTf = [R T; 0 0 0 1];
% --------------------------------------------------------------------
h = figure;hold on;grid on;
set(gcf,'Position',[0 0 1600 800], 'color','k');
set(gca,'Position',[0.01 0.01 0.99,0.99], 'color','k');
for nFrm = 10 : 10 : Len
    tic;
    % --------------------------------------------------------------------
    str = sprintf('pcd%06d.pcd', nFrm);
    DataDir = fullfile( HDLFolder, str );
    if ~exist(DataDir)
        continue;
    end
    ptOrignal = pcread(DataDir);
    % --------------------------------------------------------------------
    ptCloudIn = pcdownsample(ptOrignal,'gridAverage',0.1); % for more robust estimate
    [EffIdx, NffIdx, UnkownIdx] = GrdSegmentFun(ptCloudIn.Location', [0.40 Inf], 0);
    PtNoGrd = pointCloud(ptCloudIn.Location(NffIdx,:));
    GrdPt = ptCloudIn.Location(EffIdx,:);
    [PtCylinder,PtVehicle,PtWallPt,ValidPt] = RecoCylinderFun( PtNoGrd,CylinderFv,WallFv,VehicleFv );
    AllValidPt = pointCloud([GrdPt;ValidPt]);
    ptCloud_HDL = pcdownsample(AllValidPt,'random',0.06); % control the number of points
    % --------------------------------------------------------------------
    Tf0 = affine3d(TTf(:, :, end)');
    [tform, movingReg, rmse] = IcpFun(ptCloud_HDL, ptCloud_Map, KdMap, 30, Tf0);
    tmpTf = tform.T';
    % --------------------------------------------------------------------
    DistThr = 0.2;
    [Ratio,Rmse] = StatisticalRatioFun( KdMap, ptOrignal.Location', tmpTf, DistThr );
    % --------------------------------------------------------------------
    TRatio(end+1,:) = Ratio;
    vRmse(end+1,:) = [Rmse];
    TTf(:, :, end+1) = tmpTf;
    AssumTime(end+1,:) = toc;
    vPtNum(end+1,:) = size(ptCloud_HDL.Location, 1);
    % --------------------------------------------------------------------
    str = sprintf('Frame = %04d/%04d, PtsNum = %06d/%06d, Time = %04dms, Rato = %.4f, Rmse = %.4f', nFrm, Len...
        , size(ptCloud_HDL.Location, 1), size(ptCloud_Map.Location, 1), ceil(1000.0*AssumTime(end)),Ratio, Rmse);
    disp(str);
    % --------------------------------------------------------------------
    IS_SHOW = 1;
    if IS_SHOW
        figure(h);
        cla; hold on;  grid on;   view(3);
        TFGrdPt = myRigidTransform(GrdPt,tform);
        TFValidPt = myRigidTransform(ValidPt,tform);
        pcshow(FixPt,'k');
        pcshow(TFGrdPt,'r');
        pcshow(TFValidPt,'r');
        pcshow(ptCloudB_Reduced);
        Traj = squeeze(TTf(:, 4, 1:end));
        plot3(Traj(1, 1), Traj(2, 1), Traj(3, 1), 'rp', 'markersize', 7 );
        plot3(Traj(1, :), Traj(2, :), Traj(3, :), 'yo-', 'markersize', 5 );
        plot3(Traj(1, end), Traj(2, end), Traj(3, end), 'bp', 'markersize', 7 );
        str = sprintf('id = %03d, ratio = %.2f', nFrm, TRatio(end) );
        title(str);
        pause(0.05);
    end
    % --------------------------------------------------------------------
end
% --------------------------------------------------------------------
figure;hold on;grid on;
plot(TRatio,'r.-');
title('Ratio');
figure;hold on;grid on;
plot(vRmse,'b.-');
title('All rmse');
figure;hold on;grid on;
plot(vPtNum,'g.-');
title('Number of points');
% --------------------------------------------------------------------
TTfToEulAngShow( TTf, 1 );
% --------------------------------------------------------------------