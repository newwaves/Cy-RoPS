function B = myRigidTransform(A, TransFR, TransFT)
% A : input N X 3
% affine3d input  ([R,zeros(3,1);T',1])
if nargin == 2
    R = TransFR.T(1:3, 1:3)';
    T = TransFR.T(4, 1:3)';
end
% <R,T> input R£º 3X3  T: 3X1  ([R,T;0 0 0 1])
if nargin == 3
    R = TransFR;
    T = TransFT;
end
    B = A * R';
    B(:,1) = B(:,1) + T(1);
    B(:,2) = B(:,2) + T(2);
    B(:,3) = B(:,3) + T(3);
end