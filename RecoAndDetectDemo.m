%-------------------------------------------------------------------------------
%   Recognize the cylinder with Cy-RoPS;
%   1-Delete ground --> 2-Feature Extraction --> 3-Delete Obstacles
%   modify by Z.X.Tao E_mail: newwaves@stu.xjtu.edu.cn-- 2018-02-01 
%-------------------------------------------------------------------------------
clc; close all;clear all;
isShow = 0;
HDLFolder = '.\Data';
Len = length(dir(HDLFolder)) * 10;
addpath('.\Cy-RoPS');
%-------------------------------------------------------------------------------
load 'CylinderFv.mat';
CylinderFv = mean(cc);
load 'WallFv.mat';
WallFv = mean(cc);
load 'VehicleFv.mat';
VehicleFv = mean(cc);
%-------------------------------------------------------------------------------
h = figure;
for nFrm = 10 : 10 : Len
%-------------------------------------------------------------------------------
    str = sprintf('pcd%06d.pcd', nFrm);
    DataDir = fullfile( HDLFolder, str );
    if ~exist(DataDir)
        continue;
    end
    ptCloudIn = pcread(DataDir);
%-------------------------------------------------------------------------------
    [EffIdx, NffIdx, UnkownIdx] = GrdSegmentFun(ptCloudIn.Location', [0.40 Inf], 0);
    DelGround = ptCloudIn.Location(NffIdx,:);
%-------------------------------------------------------------------------------
    PtNoGrd = pointCloud(DelGround);
    PtNoGrd = pcdownsample(PtNoGrd,'gridAverage',0.1);
%-------------------------------------------------------------------------------
	[PtCylinder,PtVehicle,PtWallPt,PtValid] = RecoCylinderFun( PtNoGrd,CylinderFv,WallFv,VehicleFv );
%-------------------------------------------------------------------------------
    cla(h);
    hold on;grid on;
    pcshow(PtNoGrd.Location(PtCylinder,:),'g');
    pcshow(PtNoGrd.Location(PtWallPt,:),'g');
    str = sprintf('Detection cylinder - Frame = %04d',nFrm);
    title(str);
    pause(0.1);
end
%-------------------------------------------------------------------------------
