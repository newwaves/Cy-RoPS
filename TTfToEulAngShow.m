function TTfToEulAngShow( TTf, isShow )
%   analysis the euler and save the euler angles
EulerAngle = [];
DegAngle = [];
for i = 1 : 1 : length(TTf)
    pose = TTf(:, :, i);
    R = pose(1:3, 1:3);
    EulerAngle(end+1,:) = rotm2eul(R);
    DegAngle(end+1,:) = rad2deg(EulerAngle(end, :));
end
if isShow
    figure;
    hold on;
    grid on;
    plot(EulerAngle(:,1),'r');
    plot(EulerAngle(:,2),'g');
    plot(EulerAngle(:,3),'b');
    legend('Z-axis','Y-axis','X-axis');
    figure;
    hold on;
    grid on;
    plot(DegAngle(:,1),'r');
    plot(DegAngle(:,2),'g');
    plot(DegAngle(:,3),'b');
    legend('Z-axis','Y-axis','X-axis');
end

end

