function [subRoPS,GridIm] = SubSC_RoPSFun(PointsXY,XBinSize,YBinSize,isZdirection)
% this code written according to the follow paper:
% Y. Guo, F. A. Sohel, M. Bennamoun, M. Lu, and J. Wan, ¡°Rotational
% projection statistics for 3d local surface description and object recognition,¡±
% International Journal of Computer Vision, vol. 105, no. 1, pp. 63¨C86, 2013.

if nargin < 4
    isZdirection = 0;
end
NumPt = length(PointsXY);
ProjImage = zeros(XBinSize,YBinSize);
% ---------------------------------------------------
if isZdirection == 0
    minX = min(PointsXY(:,1));
    stepX = (max(PointsXY(:,1))-minX)/XBinSize;
    minY = min(PointsXY(:,2));
    stepY = (max(PointsXY(:,2))-minY)/YBinSize;
else
    minX = min(PointsXY(:,1));
    stepX = (max(PointsXY(:,1)) - minX)/XBinSize;
    minY = min(PointsXY(:,2));
    maxY = max(PointsXY(:,2));
    cenY = mean(PointsXY(:,2));
    H = 2.5;    % 2.5m
    if (maxY - minY)  < 2.5      
        Delta = ( H - (maxY - minY) ) / 2;
        maxY = maxY + Delta;
        minY = minY - Delta;
    end
    if (maxY - minY)  > 2.5   
        Delta = ( (maxY - minY) - H ) / 2;
        maxY = maxY - Delta;
        minY = minY + Delta;
    end
    stepY = (maxY - minY)/YBinSize;
end
% ---------------------------------------------------
if stepX==0 || stepY==0
    subRoPS = [0,0,0,0,0];
    return;
end

for k=1:NumPt
    IndX = ceil((PointsXY(k,1) - minX)/stepX);
    IndY = ceil((PointsXY(k,2) - minY)/stepY);
    if IndX>XBinSize      IndX = XBinSize;  end
    if IndX<1                IndX = 1;            end
    if IndY>YBinSize     IndY = YBinSize;   end
    if IndY<1                IndY = 1;            end
    ProjImage(IndX,IndY) = ProjImage(IndX,IndY)+1;
end
ProjImage = ProjImage/NumPt;

% calculate the moment of project image
meanX = 0;
meanY = 0;
entropy = 0;
for IndX = 1:XBinSize
    for IndY = 1:YBinSize
        meanX = meanX+IndX*ProjImage(IndX,IndY);
        meanY = meanY+IndY*ProjImage(IndX,IndY);
        if ProjImage(IndX,IndY)>0
            entropy = entropy - ProjImage(IndX,IndY)*log2(ProjImage(IndX,IndY));
        end
    end
end
u11 = 0;
u21 = 0;
u12 = 0;
u22 = 0;

for IndX = 1:XBinSize
    for IndY = 1:YBinSize
        u11 = u11+(IndX-meanX)*(IndY-meanY)*ProjImage(IndX,IndY);
        u21 = u21+(IndX-meanX)^2*(IndY-meanY)*ProjImage(IndX,IndY);
        u12 = u12+(IndX-meanX)*(IndY-meanY)^2*ProjImage(IndX,IndY);
        u22 = u22+(IndX-meanX)^2*(IndY-meanY)^2*ProjImage(IndX,IndY);
    end
end
subRoPS = [u11,u21,u12,u22,entropy];
if nargout > 1
    GridIm = ProjImage;
end
end