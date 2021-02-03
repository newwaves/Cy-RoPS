function [tform, movingReg, rmse] = IcpFun(ptCloudA, ptCloudB, KDMap, MaxIter, InitTF, TrimDist )
    if nargin < 5
        InitTF = affine3d();
        TrimDist = 100;
    end
    if nargin < 6
        TrimDist = 100;
    end
    if ~(isRigidTransform(InitTF))
        error(message('vision:pointcloud:rigidTransformOnly'));
    end

    tolerance = [0.01, 0.009];
    %%
    Rs = zeros(3, 3, MaxIter+1);
    Ts = zeros(3, MaxIter+1);
    % Quaternion and translation vector
    qs = [ones(1, MaxIter+1); zeros(6, MaxIter+1)];
    % RMSE
    Err = zeros(MaxIter+1, 1); 
    Rs(:,:,1) = InitTF.T(1:3, 1:3)';
    Ts(:,1) = InitTF.T(4, 1:3)';
    qs(:,1) = [rotm2quat(Rs(:,:,1))'; Ts(:,1)];
    locA = ptCloudA.Location;
    if qs(1) ~= 0 || any(qs(2:end,1))
        locA = rigidTransform(ptCloudA.Location, Rs(:,:,1), Ts(:,1));
    end
    StopIter = MaxIter;
  
    % Start ICP iterations
    for i = 1 : MaxIter
        % Find the correspondence
        [indices, dists] = knnsearch(KDMap, locA, 'dist','euclidean');
        % Remove outliers      
        keepInlierA = false(ptCloudA.Count, 1); 
        [~, idx] = sort(dists);
        upperBound = length(find(dists < TrimDist));
        keepInlierA(idx(1:upperBound)) = true;       
        IdxA = find(keepInlierA);
        IdxB = indices(keepInlierA);
        inlierDist = dists(keepInlierA);

        if numel(IdxA) < 3
            error(message('vision:pointcloud:notEnoughPoints'));
        end
        if i == 1
            Err(i) = sqrt(sum(inlierDist)/length(inlierDist));
        end

        [R, T] = pointToPlaneMetric1008(locA(IdxA, :),ptCloudB.Location(IdxB, :), ptCloudB.Normal(IdxB, :));
        
        if any(isnan(T))||any(isnan(R(:)))
            error(message('vision:pointcloud:singularMatrix'));
        end
        
        % Update the total transformation
        Rs(:,:,i+1) = R * Rs(:,:,i);
        Ts(:,i+1) = R * Ts(:,i) + T;
        
        % RMSE
        locA = rigidTransform(ptCloudA.Location, Rs(:,:,i+1), Ts(:,i+1));
        SquError = sum((locA(IdxA, :) - ptCloudB.Location(IdxB, :)).^2, 2);
        Err(i+1) = sqrt(sum(SquError)/length(SquError));

        qs(:,i+1) = [rotm2quat(Rs(:,:,i+1))'; Ts(:,i+1)];
        % Compute the mean difference in R/T from the recent three iterations.
        [dR, dT] = getChangesInTransformation;
        % Stop ICP if it already converges
        if dT <= tolerance(1) && dR <= tolerance(2)
            StopIter = i;
            break;
        end
    end
    % Make the R to be orthogonal as much as possible
    R = Rs(:,:,StopIter+1)';
    [U, ~, V] = svd(R);
    R = U * V';
    tformMatrix = [R, zeros(3,1);Ts(:, StopIter+1)',  1];
    tform = affine3d(tformMatrix);
    rmse = Err(StopIter+1);
    if nargout >= 2
        movingReg = pctransform(ptCloudA, tform);
    end
    %======================================================================
    function [dR, dT] = getChangesInTransformation
        dR = 0;
        dT = 0;
        count = 0;
        for k = max(i-2,1):i
            % Rotation difference in radians
            rdiff = acos(dot(qs(1:4,k),qs(1:4,k+1))/(norm(qs(1:4,k))*norm(qs(1:4,k+1))));
            % Euclidean difference
            tdiff = sqrt(sum((Ts(:,k)-Ts(:,k+1)).^2));
            dR = dR + rdiff;
            dT = dT + tdiff;
            count = count + 1;
        end
        dT = dT/count;
        dR = dR/count;
    end
end

%==========================================================================
% Determine if transformation is rigid transformation
%==========================================================================
function tf = isRigidTransform(tform)
    singularValues = svd(tform.T(1:tform.Dimensionality,1:tform.Dimensionality));
    tf = max(singularValues)-min(singularValues) < 100*eps(max(singularValues(:)));
    tf = tf && abs(det(tform.T)-1) < 100*eps(class(tform.T));
end 
%%
%==========================================================================
function B = rigidTransform(A, R, T)
    B = A * R';
    B(:,1) = B(:,1) + T(1);
    B(:,2) = B(:,2) + T(2);
    B(:,3) = B(:,3) + T(3);
end
%%
function [R, T] = pointToPlaneMetric1008(p, q, nv)
    % Set up the linear system
    cn = [cross(p,nv,2),nv];
    C = cn'*cn;
    qp = q-p;
    b =  [sum(sum(qp.*repmat(cn(:,1),1,3).*nv, 2));
          sum(sum(qp.*repmat(cn(:,2),1,3).*nv, 2));
          sum(sum(qp.*repmat(cn(:,3),1,3).*nv, 2));
          sum(sum(qp.*repmat(cn(:,4),1,3).*nv, 2));
          sum(sum(qp.*repmat(cn(:,5),1,3).*nv, 2));
          sum(sum(qp.*repmat(cn(:,6),1,3).*nv, 2))];
    % X is [alpha, beta, gamma, Tx, Ty, Tz]
    
    X = C\b;  % 可能 矩阵接近奇异值
   
    cx = cos(X(1)); 
    cy = cos(X(2)); 
    cz = cos(X(3)); 
    sx = sin(X(1)); 
    sy = sin(X(2)); 
    sz = sin(X(3)); 
    %  euler to rotation matrix -> [R,T]
    R = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
         cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
           -sy,          sx*cy,          cx*cy];

    T = X(4:6);
end