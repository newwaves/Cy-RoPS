function [SimpleRoPS,FeatHist] = CalSC_RoPS(pcSeed,PtData,Radius,BinSize,KDpcData)
% this code written according to the follow paper:
% Y. Guo, F. A. Sohel, M. Bennamoun, M. Lu, and J. Wan, ¡°Rotational
% projection statistics for 3d local surface description and object recognition,¡±
% International Journal of Computer Vision, vol. 105, no. 1, pp. 63¨C86, 2013.

    if nargin < 5
        KDpcData = KDTreeSearcher(PtData(:,1:2),'Distance','euclidean');
    end
    if nargin < 4
        BinSize = 10;
        KDpcData = KDTreeSearcher(PtData(:,1:2),'Distance','euclidean');
    end
    if nargin < 3
        Radius = 0.5;
        BinSize = 10;
        KDpcData = KDTreeSearcher(PtData(:,1:2),'Distance','euclidean');
    end

    [neighborIdx2,~]= rangesearch(KDpcData,pcSeed(1:2),Radius);    
    neighborIdx2 = cell2mat(neighborIdx2);
    neighborIdx = neighborIdx2(2:end); % delete self-data
    neighbNum = length(neighborIdx);
    
    rotaSize = 3;
    interval = pi/2/rotaSize;
    
    if neighbNum<=1
        SimpleRoPS = ones(1,rotaSize*15)/sum(ones(1,rotaSize*15));
        return;
    end
    isShow = 0;
    isAddLRF = 1;   % Note : 0 for FeatureSimilarity.m file
    if isAddLRF
        %  estimate principal normal
        TmpPt = PtData(neighborIdx,:);
        %--------------------------------------------------
        % construct Local Reference Frame
        %--------------------------------------------------
        myCov = cov(TmpPt);
        [U,S,V] = svd(myCov);
        AppxNorm = U(:,3)'; % along the normal direction      
        zAxis = [0,0,1];    % modifyed the principal direction
        zAxis = zAxis / norm(zAxis);
        yAxis = cross(AppxNorm,zAxis);
        yAxis = yAxis / norm(yAxis);
        xAxis = cross(yAxis,zAxis);
        xAxis = xAxis / norm(xAxis);
        if dot(([0,0,0]-pcSeed),xAxis) < 0   % mean direct to original point
            xAxis = -xAxis; % correspond to mininum eigen value
        end
        %--------------------------------------------------
        Fuv = [xAxis;yAxis;zAxis];
        Fwx = [1,0,0;0,1,0;0,0,1];
        TF = Fwx / Fuv; % means - Fwx * inv(Fuv) : coordinate transform
        neighbor = (PtData(neighborIdx,:) - pcSeed) * TF;
        if isShow
            figure;hold on;grid on;
            axis([-1,1,-1,1,-1.5,1.5]);
            pcshow(TmpPt,'k');
            ShowCylinderOnPC(Radius,pcSeed(1),pcSeed(2),1.5,-1.5);
            ShowFrameSysFun(pcSeed,Fuv);
            pcshow(neighbor,'r');
            ShowCylinderOnPC(Radius,0,0,2.5);
            ShowFrameSysFun([0,0,0],Fwx);
        end
    else
        neighbor = PtData(neighborIdx,:);
    end
    RoPS = [];

    %calculate the sub-feature of the keypoint along the Z axis
        
    for rotaIdx = 1:rotaSize
        rotaAngle = (rotaIdx-1)*interval + interval/2;
        R = [cos(rotaAngle) sin(rotaAngle) 0; -sin(rotaAngle) cos(rotaAngle) 0; 0 0 1]';
        rotaNeighbor = neighbor*R;
        %projection on the XY plane 
        projNeighborXY = [rotaNeighbor(:,1),rotaNeighbor(:,2)];
        XbinSize = 5;
        YbinSize = 5;
        histTemp = SubSC_RoPSFun(projNeighborXY,XbinSize,YbinSize);
        RoPS = [RoPS,histTemp];       
        %projection on the XZ plane 
        projNeighborXZ = [rotaNeighbor(:,1),rotaNeighbor(:,3)];
        XbinSize = 5;
        YbinSize = 10;
        histTemp = SubSC_RoPSFun(projNeighborXZ,XbinSize,YbinSize,1);
        RoPS = [RoPS,histTemp];
        %projection on the YZ plane 
        projNeighborYZ = [rotaNeighbor(:,2),rotaNeighbor(:,3)];
        XbinSize = 5;
        YbinSize = 10;
        histTemp = SubSC_RoPSFun(projNeighborYZ,XbinSize,YbinSize,1);
        RoPS = [RoPS,histTemp];
    end
    SimpleRoPS = (RoPS)/sum(RoPS);
    if nargout > 1
        [FeatHist,~] = hist(SimpleRoPS,BinSize); 
        FeatHist = FeatHist / sum(FeatHist);
    end 
end