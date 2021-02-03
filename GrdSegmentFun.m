function [EffIdx, NffIdx, UnkownIdx] = GrdSegmentFun(pcData, GapThr, IS_SHOW)
% this code from : https://github.com/DrGabor/LiDAR
    if nargin == 0
        clc; close all;
        DataFolder = 'F:\DATA3DFROMGABOR\Record-2016-10-24-10-54-01(HighWayL)\BinaryData';
        nFrm = 1000; 
        DataDir = fullfile(DataFolder, sprintf('Binary%d.txt', nFrm)); 
        pcData = HDLAnalyserNew(DataDir); 
        pcData = pcData(1:3, :);
        GapThr = [0.10 Inf];
        IS_SHOW = 1; 
    end
    % addpath('C:\Program Files\MATLAB\R2017a\bin\Ransac\'); 
    tic;
    RadArray = 0.0 : 1.0 : 80.0;
    AngRes = deg2rad(5.0);
    [SegID, BinID, EffIdx] = CvtPtsToPolar(pcData, RadArray, AngRes);
    SegID = SegID(EffIdx);
    BinID = BinID(EffIdx);
    data = pcData(:, EffIdx);
    maxAngLen = ceil(2*pi/AngRes);
    maxBinLen = length(RadArray);
    ind = sub2ind( [maxAngLen maxBinLen], SegID, BinID );
    [B I] = sort(ind);
    [C ia ic] = unique(B);
    ss = struct('Points', [], 'RawIdx', [], 'Gap', -Inf, 'MPts', -Inf(3, 1), 'CoorI', [] );
    TotalPixel = [];
    [seg bin] = ind2sub([maxAngLen maxBinLen], C);
    for i = 1 : 1 : length(C)   % C = ind(ia) or
        id0 = ia(i);
        if i == length(C)
            id1 = length(ind);
        else
            id1 = ia(i+1) - 1;
        end
        Idx = I(id0:1:id1);
        if length(Idx) > 10
            bTest = 1;
        end
        a = unique(ind(Idx));
        if length(a) ~= 1
            error('wrong!');
        end
        tmp = ss;
        tmp.RawIdx = EffIdx(Idx);
        tmp.CoorI = [seg(i); bin(i)];
        TotalPixel = [TotalPixel tmp];
    end
    CoorI = cat(2, TotalPixel(:).CoorI);
    ObsIdx = [];
    UnkownSeg0 = [];
    ss = struct('data', [] );
    RMInfo = repmat(ss, maxAngLen, maxBinLen);
    UnkownIdx = [];
    for AngID = 1 : 1 : maxAngLen
        if AngID == 15
            bTest = 1; 
        end
        SelIdx = find(CoorI(1, :) == AngID );
        SelPts = [];
        for i = 1 : 1 : length(SelIdx)
            tmp = TotalPixel(SelIdx(i));
            Pts = pcData(:, tmp.RawIdx);
            [~, id] = min(Pts(3, :));
            SelPts(:, end+1) = Pts(:, id);
            bTest = 1;
        end
        EIdx = cat(2, TotalPixel(SelIdx).RawIdx);
        if size(SelPts, 2) < 10
            UnkownSeg0 = [UnkownSeg0 AngID];
            str = sprintf('AngID = %03d/%03d, too few points, continue!', AngID, maxAngLen);
            continue;
        end
        %%
        Dist = sqrt( SelPts(1, :).^2 + SelPts(2, :).^2 );
        Pts = [Dist; SelPts(3, :)];
        
        [V, L, inliers] = ransacfit2Dline(Pts, 0.1, 0);
        Dist = ( Pts(1, :)*V(1) + Pts(2, :)*V(2)+V(3) ) / norm(V(1:2));
        %%%%% y = kx + b.
        k = -V(1)/V(2);
        b = -V(3)/V(2);

        H0 = k*5.0+b;
        idx = find(abs(Dist) <= 0.10);

        if length(idx) < 10 || ( H0 > -1.5 || H0 < -2.1)  % the orignal ground height is -1.8m.
            UnkownSeg0 = [UnkownSeg0 AngID];
            if length(idx) < 10
                str = sprintf('AngID = %03d/%03d, too few initial points, continue!', AngID, maxAngLen);
            else
                str = sprintf('AngID = %03d/%03d, ransac finds a wrong line, continue!', AngID, maxAngLen);
            end
            continue;
        end
        str = sprintf('AngID = %03d/%03d', AngID, maxAngLen);

        tData = 2.0;
        tDist = 0.15;

        Hyp.cov = log([6.0 sqrt(1.3298)]);
        Hyp.mean = [0 0];
        Hyp.lik = log(0.1);
        out = GPRGrdSegFun_WithFixedHyp(Pts, idx, tData, tDist, Hyp); 
        % out = GPRGrdSegFun(Pts, idx', tData, tDist, [], 0);
        VIdx = out.ValidIdx;
        xs = Pts(1, :);
        ys = Pts(2, :);
        data_t = out.data(out.ValidIdx, :); 
        [mu, s2, F, EffIdx, NffIdx] = GPR_PredictFun(out.hyp, data_t, xs, ys, tData, tDist);

        CorrectH = ys';
        CorrectH(NffIdx) = mu(NffIdx);
        % CorrectH = mu;
        for i = 1 : 1 : length(SelIdx)
            tmp = TotalPixel(SelIdx(i));
            Pts = pcData(:, tmp.RawIdx);
            tmpDiff = Pts(3, :) - CorrectH(i);
            idx = find(tmpDiff >= GapThr(1) & tmpDiff <= GapThr(2) );
            ObsIdx = [ObsIdx tmp.RawIdx(idx)];
            tmpIdx = find(tmpDiff > GapThr(2));
            UnkownIdx = [UnkownIdx tmp.RawIdx(tmpIdx)];
            bTest = 1;
            if size(Pts, 2) > 1
                MPt = mean(Pts')';
            else
                MPt = Pts;
            end
            MPt(3) = CorrectH(i);
            RMInfo( tmp.CoorI(1), tmp.CoorI(2) ).data = MPt;
        end
        bTest = 1;
    end
    %%%%%%%%%%% check the unkown segments.  
    GrdPts = cat(2, RMInfo(:).data);
    GrdMd = createns(GrdPts');
    UnkownSeg = [];

    NewObsIdx = [];
    rArray = [2.0 5.0 10.0 15.0 20.0 25.0 30.0];
    for id = 1:1:length(UnkownSeg0)
        AngID = UnkownSeg0(id);
        SelIdx = find(CoorI(1, :) == AngID );
        TTIdx = cat(2, TotalPixel(SelIdx).RawIdx);
        SelPts = [];
        for i = 1 : 1 : length(SelIdx)
            tmp = TotalPixel(SelIdx(i));
            Pts = pcData(:, tmp.RawIdx);
            if size(Pts, 2) > 1
                MPt = mean(Pts')';
            else
                MPt = Pts;
            end
            for rId = 1:1:length(rArray)
                [NNIdx, DD] = rangesearch(GrdMd, MPt', rArray(rId));
                if ~isempty(NNIdx{1})
                    break;
                end
            end
            if isempty(NNIdx{1})
                UnkownIdx = [UnkownIdx tmp.RawIdx];
                UnkownSeg = [UnkownSeg AngID];
                continue;
            end
            NNIdx = NNIdx{1};
            tmpPts = GrdPts(:, NNIdx);
            H = mean(tmpPts(3, :));
            tmpDiff = Pts(3, :) - H;
            idx = find(tmpDiff >= 0.3 & tmpDiff <= GapThr(2) );
            NewObsIdx = [NewObsIdx tmp.RawIdx(idx)];
            tmpIdx = find(tmpDiff > GapThr(2));
            UnkownIdx = [UnkownIdx tmp.RawIdx(tmpIdx)];
        end
        bTest = 1;
    end
    ObsIdx = [ObsIdx NewObsIdx];    
    Dist = sqrt(pcData(1, :).^2 + pcData(2, :).^2);
    NIdx = find(Dist > RadArray(end) );
    GrdIdx = find(~ismember(1:1:size(pcData,2), ObsIdx) & ~ismember(1:1:size(pcData,2), UnkownIdx) & ~ismember(1:1:size(pcData,2), NIdx) );
    
    if IS_SHOW
        figure;
        hold on;
        grid on;
        view(3);
        if ~isempty(GrdIdx)
            showPointCloud(pcData(1:3, GrdIdx)', 'g');
        end
        if ~isempty(UnkownIdx)
            showPointCloud(pcData(1:3, UnkownIdx)', 'k');
        end
        if ~isempty(ObsIdx)
            showPointCloud(pcData(1:3, ObsIdx)', 'r');
        end
        title('Segmentation Results');
    end
    toc
    NffIdx = ObsIdx;
    EffIdx = GrdIdx;
    bTest = 1;
end

%% --------------------------------------------------------------------
function [varargout]  = CvtPtsToPolar( pcData, varargin )
    switch (nargin-1)
        case 0
            RadArray = [ 0.0 : 0.2 : 20.0 20.5 : 0.5 : 50.0 ];
            AngRes   = deg2rad(1.0);
        case 1
            AngRes = varargin{1};
            RadArray = [ 0.0 : 0.2 : 20.0 20.5 : 0.5 : 50.0 ];
        case 2
            RadArray = varargin{1};
            AngRes   = varargin{2};
        otherwise
            error('Invalid input!\n');
    end
    Angle = wrapTo2PiFun( atan2( pcData(2, :) , pcData(1, :) ) );
    SegID = ceil( (Angle+1e-6)/ AngRes );
    maxLen = ceil(2*pi/AngRes); 
    SegID(SegID>maxLen) = maxLen;
    USE_ANGLE_RANGE = 0; 
    if USE_ANGLE_RANGE
        EffAngIdx = find( Angle <= deg2rad(90.0) | Angle >= deg2rad(270.0) ); 
    end
    Radius = sqrt( pcData(1, :).^2 + pcData(2, :).^2 );
    BinID = -1 * ones(1, size(pcData, 2) );   % BinID = 1 stores invalid points, such as radius smaller than 4.5m or bigger than 5.0m.
    for i = 1 : 1 : (length(RadArray)-1)
        Idx = find( Radius >= RadArray(i) & Radius < RadArray(i+1) );
        if USE_ANGLE_RANGE
            ind = ismember(Idx, EffAngIdx); 
            Idx = Idx(ind); 
        end
        BinID(Idx) = i;
    end
    if nargout == 2
        varargout{1} = SegID;
        varargout{2} = BinID;
    end
    if nargout == 3
        varargout{1} = SegID;
        varargout{2} = BinID;
        varargout{3} = find(BinID ~= -1);
    end
end
%% --------------------------------------------------------------------
function [ AngleNew ] = wrapTo2PiFun(Angle)
AngleNew = Angle; 
idx = find(Angle < 0 );
AngleNew(idx) = AngleNew(idx) + 2*pi; 
end
%% --------------------------------------------------------------------
function [ymu, yms2] = myGPRFun(hyp, x, y, xs)
    IS_SHOW = 0; 
    if nargin == 0
        hyp = [];
        hyp.mean = [0.1 0.2];
        hyp.cov = log([1.2 1.5]);
        hyp.lik = log(0.1);
        x = [-5:0.1:5]';
        [Kff, ~] = CalKFun(x, x, hyp.cov);
        Mx = hyp.mean(1)*x + hyp.mean(2); 
        y = mvnrnd(Mx, Kff); 
        y = y'; 
        xs = [-10:0.1:10]'; % 10*rand(50, 1)-5;
        IS_SHOW = 1; 
    end
    m0 = hyp.mean(1);
    m1 = hyp.mean(2);
    sn = exp(hyp.lik);
    sn2 = sn*sn;
    n = length(x);
    p = hyp.cov;
    [Kuu, ~] = CalKFun(xs, xs, p);
    [Kuf, ~] = CalKFun(xs, x, p);
    [Kff, ~] = CalKFun(x, x, p);
    A = Kff + sn2*eye(n); 
    L = chol(A, 'lower');
    iL = inv(L);
    iA = iL'*iL;
    % tt = diag(iA*A); 
    % iA = inv(A); 
    ymu = m0*xs + m1 + Kuf*iA*(y - m0*x-m1);
    m = length(xs);
    fs2 = diag(Kuu - Kuf*iA*Kuf');
    yms2 = fs2+sn2;
    bTest = 1;
    if IS_SHOW
        Mu = ymu;
        S2 = yms2;
        Xs = xs;
        tData = 2.0; 
        a = Mu+tData*sqrt(S2); 
        b = flipdim(Mu-tData*sqrt(S2),1); 
        F = [a;b];
        figure;
        hold on;
        grid on;
        axis equal;
        %%%%%%%%%% draw envelop.
    %     plot(xs, a, 'r.--'); 
    %     plot(xs, b, 'g.--'); 
        fill([Xs; flipdim(Xs,1)], F, [7 7 7]/8)
        plot(Xs, Mu, 'b.');
        plot(x, y, 'r.'); 
    end
end
%% --------------------------------------------------------------------
function [K, dK] = CalKFun(x, y, p)
l = exp(p(1));
sf = exp(p(2));
l2 = l*l;
sf2 = sf*sf;
m = length(x);
n = length(y);
dK = {};
X = repmat(x, 1, n); 
Y = repmat(y', m, 1); 
R = (X - Y) .*(X - Y);  
K = sf2*exp(-0.5/l2*R); 
dK{1} = 1/l2 * R.*K;
dK{2} = 2*K;
end
%% --------------------------------------------------------------------
function [mu, s2, F, EffIdx, NffIdx] = GPR_PredictFun(hyp, data, xs, ys, tData, tDist)
    if size(data, 2) ~= 2
        data = data';
    end
    if size(data, 2) ~= 2
        error('Input for GPR_PredictFun is wrong!');
    end
    if isrow(xs)
        xs = xs'; 
    end
    if isrow(ys)
        ys = ys'; 
    end
    x = data(:, 1);
    y = data(:, 2);
    [mu s2] = myGPRFun(hyp, x, y, xs);
    F = [mu+tData*sqrt(s2); flipdim(mu-tData*sqrt(s2),1)];
    if ~isempty(ys)
        Ratio = (mu - ys)./sqrt(s2) / tData;
        %%%%%%%%% update idx.
        EffIdx = find( abs(Ratio) <= 1 & sqrt(s2) <= tDist );
        NffIdx = find( ~ismember( 1:1:length(xs), EffIdx) );
        if isrow(NffIdx)
            NffIdx = NffIdx';
        end
    else
        EffIdx = [];
        NffIdx = [];
    end
end
%% --------------------------------------------------------------------
function out = GPRGrdSegFun_WithFixedHyp(data, VIdx, tData, tDist, Hyp)
    if nargin == 0
    end
    if size(data, 2) ~= 2
        data = data';
    end
    if size(data, 2) ~= 2
        error('GPRGrdSegFun_WithFixedHyp() is wrong!');
    end
    out = []; 
    out = struct('data', [], 'ValidIdx', [], 'NIdx', [], 'hyp', [] );
    out.data = data;
    nLen = max(size(data)); 
    NIdx = find( ~ismember([1:1:nLen], VIdx) ); 
    for i = 1 : 1 : length(NIdx)
        id = NIdx(i);
        pt = data(id, :);
        xt = data(VIdx, 1); 
        yt = data(VIdx, 2); 
        [mu s2] = myGPRFun(Hyp, xt, yt, pt(1));
        Ratio = (mu - pt(2))./sqrt(s2) / tData;
        %%%%%%%%% update idx.
       if abs(Ratio) <= 1 & sqrt(s2) <= tDist 
           VIdx(end+1) = id; 
       end
    end
    NIdx = find( ~ismember([1:1:nLen], VIdx) ); 
    out.ValidIdx = VIdx;
    out.NIdx = NIdx;
    out.hyp  = Hyp;
end
%% --------------------------------------------------------------------