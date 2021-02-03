function [Ratio,Rmse] = StatisticalRatioFun( Md, MovData, Tf, DistThreshold )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% Find the correspondence
    Md_Ref = Md; % createns(RefData');
    R = Tf(1:3, 1:3); 
    T = Tf(1:3, end); 
    
    %TransF = affine3d([R,zeros(3,1);T',1]);
    %HdlData = pctransform(HdlData,TransF);
    AftData = Loc2Glo(MovData, R', T);
    [~, dists] = knnsearch(Md_Ref, AftData');
    ValidIdx = find(dists <= DistThreshold);
    nInlier = length(ValidIdx);
    Ratio = nInlier/size(MovData, 2);
    if nargout > 1
        Rmse = mean(dists(ValidIdx,:));
    end
end

