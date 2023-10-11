% plotExFreeWalkTraj2.m
%
% Script to plot example free walking trajectories (in XY plane) and
%  corresponding forward and yaw velocities
%
% CREATED: 10/11/23 - HHY
%
% UPDATED:
%   10/11/23 - HHY


figure;
plot(bodytraj.tZeroed, bodytraj.angVelSmoS);

figure;
plot(bodytraj.tZeroed, bodytraj.fwdVelSmoS);

figure;
plot(bodytraj.angVelSmoS);

figure;
plot(bodytraj.fwdVelSmoS);

% these are the current example traces
% trial-1-20181105-115917.mat
% selInd = 40207:40276;
% selInd = 37434:37502;
% selInd = 41025:41115;
% selInd = 40698:40770;
% selInd = 40216:40275;
% selInd = 39637:39698;
% selInd = 41120:41182;
% selInd = 13705:13800;
% selInd = 14546:14610;
% selInd = 4122:4181;
% selInd = 1087:1141;
% selInd = 780:838;
% selInd = 947:1026;

% example trace from
load('trial-23-20181107-192031.mat')
selInd = 3080:3153;



% interpolate to 100 Hz
interpFramerate = 100;
ifi = 1/interpFramerate;

newT = bodytraj.tZeroed(1):ifi:bodytraj.tZeroed(end);

newX = interp1(bodytraj.tZeroed,bodytraj.x, newT,'linear');
newY = interp1(bodytraj.tZeroed,bodytraj.y, newT, 'linear');
newAngVel = interp1(bodytraj.tZeroed,bodytraj.angVelSmoM, newT, 'linear');
newFwdVel = interp1(bodytraj.tZeroed,bodytraj.fwdVelSmoM, newT, 'linear');
newUnWrAng = interp1(bodytraj.tZeroed,bodytraj.unWrAngSmoM, newT, 'linear');
newWrAng = wrapTo360(newUnWrAng);

% selInd = 32740:32769; % new decrease in forward velocity
% 
% selInd = 32499:32528; % increase in forward velocity with turn
% selInd = 33120:33149; % decrease in forward velocity with turn

% interpSelInd = 32456:32482;
% interpSelInd = 30215:30242;

interpSelInd = find(newT >= bodytraj.tZeroed(selInd(1)) & ...
    newT <= bodytraj.tZeroed(selInd(end)));

totT = newT(interpSelInd(end)) - newT(interpSelInd(1));


% plot angular velocity
% figure;
% plot(newT(interpSelInd), newAngVel(interpSelInd));
% xlabel('Time (s)');
% ylabel('Yaw velocity (deg/s)');
% xlim([newT(interpSelInd(1)) newT(interpSelInd(end))]);
% ylim([-50 250]);
% 
% % plot forward velocity
% figure;
% plot(newT(interpSelInd), newFwdVel(interpSelInd) * 1000);
% xlabel('Time (s)');
% ylabel('Forward velocity (mm/s)');
% xlim([newT(interpSelInd(1)) newT(interpSelInd(end))]);
% ylim([5 20]);
% % ylim([10 25]);

% compute vectors for quiver plot, angle is heading direction and length of
%  arrow is proportional to angular velocity
% theseAngs = newWrAng(interpSelInd);
% theseTan = tand(theseAngs);
% 
% theseAngVel = newAngVel(interpSelInd);
% 
% a = sqrt(theseAngVel.^2 ./ (theseTan.^2 + 1));
% o = a .* theseTan;


% plot trajectory, with arrows for heading direction and angular velocity
c = colormap('lines');
figure;
colormap('parula');
% quiver(newX(selInd)*1000,newY(selInd)*1000,a, ...
%     o);
% hold on;
scatter(newX(interpSelInd)*1000,newY(interpSelInd)*1000, [], ...
    linspace(0,totT,length(interpSelInd)), 'filled');
colorbar
axis('equal');
hold on;

% plot foot falls
% get absolute positions for when foot touches ground
minLegInds = legSteps.minIndsAll(legSteps.minIndsAll >= selInd(1) & ...
    legSteps.minIndsAll <= selInd(end));
minLegIndsWhichLeg = legSteps.minWhichLeg(legSteps.minIndsAll >= selInd(1) & ...
    legSteps.minIndsAll <= selInd(end));

for i = 1:6
    theseInd = minLegInds(minLegIndsWhichLeg == i);
    plot(legTrack.aLegX(theseInd,i)*1000, legTrack.aLegY(theseInd,i)*1000,...
        'Marker','o', 'Color', c(i,:));
    hold on;
end
axis('equal');

% plot angular velocity
figure;

subplot(3,1,1);
plot(newT(interpSelInd), newAngVel(interpSelInd));
xlabel('Time (s)');
ylabel('Yaw velocity (deg/s)');
xlim([newT(interpSelInd(1)) newT(interpSelInd(end))]);
% ylim([-50 250]);

% plot forward velocity
subplot(3,1,2);
plot(newT(interpSelInd), newFwdVel(interpSelInd) * 1000);
xlabel('Time (s)');
ylabel('Forward velocity (mm/s)');
xlim([newT(interpSelInd(1)) newT(interpSelInd(end))]);
% ylim([5 20]);
% ylim([10 25]);

% plot steps body coordinates, x
subplot(3,1,3);
plot(legTrack.t(selInd), legTrack.srnfLegX(selInd,1:6));
xlim([newT(interpSelInd(1)) newT(interpSelInd(end))]);
set(gca, 'ydir', 'reverse')
