%%
function [vCylinderInx,vVehicleInx,vWallInx, ValidPt] = RecoCylinderFun( PtData,CylinderFv, WallFv, VehicleFv )
%-------------------------------------------------------------------------------
Radius = 0.5;
Method = 'Norm';'Exp';'Cos';
switch Method
    case 'Cos'
        ThrVal = 0.9;
    case 'Exp'
        ThrVal = 0.9;
    case 'Norm'
        ThrVal = 0.1;
    otherwise
        ThrVal = 0.1;
end

PtSeed = pcdownsample( PtData,'gridAverage',1.0);
KdPtData = KDTreeSearcher(PtData.Location(:,1:2),'Distance','euclidean');
[vAllIdx,~]= rangesearch(KdPtData,PtSeed.Location(:,1:2),Radius);
vCyDist = [];
vVeDist = [];
vWaDist = [];
vValidPtIdx = {};
for i = 1 : 1 : length(vAllIdx)
    nIdx = vAllIdx{i,:};
    if length(nIdx) < 50
        continue;
    end    
    vValidPtIdx(end+1,:) = {nIdx};
    TmpPt = PtData.Location(nIdx,:);
    CtSeed = PtSeed.Location(i,:);
    SimpleFV = CalSC_RoPS(CtSeed,TmpPt,Radius);
    %------------------------------------------
    switch Method
        case 'Cos'
            vCyDist(end+1,:) = SimilarityFunCos(CylinderFv,SimpleFV);
            vVeDist(end+1,:) = SimilarityFunCos(VehicleFv,SimpleFV);
            vWaDist(end+1,:) = SimilarityFunCos(WallFv,SimpleFV);
        case 'Exp'
            vCyDist(end+1,:) = SimilarityFunExp(CylinderFv,SimpleFV);
            vVeDist(end+1,:) = SimilarityFunExp(VehicleFv,SimpleFV);
            vWaDist(end+1,:) = SimilarityFunExp(WallFv,SimpleFV);
        case 'Norm'
            vCyDist(end+1,:) = norm(CylinderFv-SimpleFV);
            vVeDist(end+1,:) = norm(VehicleFv-SimpleFV);
            vWaDist(end+1,:) = norm(WallFv-SimpleFV);
    end
end
%-------------------------------------------------------------------------------
vAllProb = [vCyDist,vVeDist,vWaDist];
vCylinderInx = [];
vVehicleInx = [];
vWallInx = [];
for i = 1 : 1 : length(vAllProb)
    nIdx = vValidPtIdx{i,:}';
    if strcmp(Method,'Norm')==1
        [minV,Idx] = min(vAllProb(i,:));
        if Idx == 1 && minV < ThrVal
            vCylinderInx = [vCylinderInx;nIdx];
        end
        if Idx == 2 && minV < ThrVal
            vVehicleInx = [vVehicleInx;nIdx];
        end
        if Idx == 3 && minV < ThrVal
            vWallInx = [vWallInx;nIdx];
        end
    else
        [maxV,Idx] = max(vAllProb(i,:));
        if Idx == 1 && maxV > ThrVal
            vCylinderInx = [vCylinderInx;nIdx];
        end
        if Idx == 2 && maxV > ThrVal
            vVehicleInx = [vVehicleInx;nIdx];
        end
        if Idx == 3 && maxV > ThrVal
            vWallInx = [vWallInx;nIdx];
        end
    end
end
vCylinderInx = unique(vCylinderInx,'stable');
vVehicleInx = unique(vVehicleInx,'stable');
vWallInx = unique(vWallInx,'stable');
%-------------------------------------------------------------------------------
if nargout > 3
    vAllInd = unique([vCylinderInx;vVehicleInx;vWallInx],'stable');
    [ValidIdx,~] = setdiff(vAllInd,vVehicleInx,'stable');
    ValidPt = PtData.Location(ValidIdx,:);
end
%-------------------------------------------------------------------------------
end