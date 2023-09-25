% plotExFreeWalkTraj.m
%
% Script to plot example free walking trajectories (in XY plane) and
%  corresponding forward and yaw velocities
%
% CREATED: 9/8/23 - HHY
%
% UPDATED:
%   9/8/23 - HHY


figure;
plot(bodytraj.tZeroed, bodytraj.angVelSmoS);

figure;
plot(bodytraj.tZeroed, bodytraj.fwdVelSmoS);

% these are the current example traces
% trial-1-20181105-115917.mat
selInd = 32499:32528; % increase in forward velocity with turn

% trial-20-20181107-191037.mat
selInd = 8040:8069; % slows down

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





% plot angular velocity
figure;
plot(newT(selInd), newAngVel(selInd));
xlabel('Time (s)');
ylabel('Yaw velocity (deg/s)');
xlim([newT(selInd(1)) newT(selInd(end))]);
ylim([-50 250]);

% plot forward velocity
figure;
plot(newT(selInd), newFwdVel(selInd) * 1000);
xlabel('Time (s)');
ylabel('Forward velocity (mm/s)');
xlim([newT(selInd(1)) newT(selInd(end))]);
ylim([5 20]);
% ylim([10 25]);

% compute vectors for quiver plot, angle is heading direction and length of
%  arrow is proportional to angular velocity
theseAngs = newWrAng(selInd);
theseTan = tand(theseAngs);

theseAngVel = newAngVel(selInd);

a = sqrt(theseAngVel.^2 ./ (theseTan.^2 + 1));
o = a .* theseTan;


% plot trajectory, with arrows for heading direction and angular velocity
figure;
% quiver(newX(selInd)*1000,newY(selInd)*1000,a, ...
%     o);
% hold on;
scatter(newX(selInd)*1000,newY(selInd)*1000, [], ...
    linspace(0,0.3,length(selInd)), 'filled');
colorbar
axis('equal');
xlim([571.5 577.5])
xlim([607 613])

xlim([223 229])

% set(gca, 'xdir', 'reverse');
set(gca, 'ydir', 'reverse')

%% not interpolated

% selInd = 40570:40617;

% selInd = 39703:39765;
% selInd = 39709:39750;
selInd = 40285:40316; % increase in forward velocity
% selInd = 39815:39890; % decrease in forward velocity
selInd = 41050:41090; % decrease in forward velocity


% figure;
% plot(bodytraj.tZeroed(selInd), bodytraj.angVelSmoS(selInd));
% 
% figure;
% plot(bodytraj.tZeroed(selInd), bodytraj.fwdVelSmoS(selInd));

figure;
plot(bodytraj.tZeroed(selInd), bodytraj.angVelSmoM(selInd));

figure;
plot(bodytraj.tZeroed(selInd), bodytraj.fwdVelSmoM(selInd));


% figure;
% scatter(bodytraj.x(selInd),bodytraj.y(selInd), [], ...
%     linspace(1,10,length(selInd)), 'filled');
% axis('equal')

theseAngs = bodytraj.wrAngSmoM(selInd);
theseTan = tand(theseAngs);

theseAngVel = bodytraj.angVelSmoM(selInd);

a = sqrt(theseAngVel.^2 ./ (theseTan.^2 + 1));
o = a .* theseTan;


figure;
% quiver(bodytraj.x(selInd),bodytraj.y(selInd),ones(size(theseTan)), ...
%     theseTan,0.5);

quiver(bodytraj.x(selInd),bodytraj.y(selInd),a, ...
    o);
hold on;
scatter(bodytraj.x(selInd),bodytraj.y(selInd), [], ...
    linspace(1,10,length(selInd)), 'filled');
axis('equal');