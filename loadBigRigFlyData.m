% loadBigRigFlyData.m
%
% Script to load big-rig data from Luke and output fly trajectories and
%  velocities 
% Adapted from main-for-helen jupyter notebook from Luke
%
% CREATED: 6/11/19
% UPDATED: 6/24/19

%% define constants
TIME_RES = 0.01; % sampling time to interpolate fly behavior to
% smoothing of fly behavior, window for smoothdata function, 
%  method gaussian
SIGMA = 10; 
%% define data paths

exptPath = '/Users/hyang/Dropbox (Personal)/Clandinin Lab/Luke_bigrig';

expts = {'exp-20181102-165424', ...
         'exp-20181102-175232', ...
         'exp-20181103-184106', ...
         'exp-20181104-162518', ...
         'exp-20181105-115608', ...
         'exp-20181107-181316', ...
         'exp-20181108-111101', ...
         'exp-20181108-143044', ...
         'exp-20181109-084314'};

curDir = pwd;
     
%% load data

% trial counter
trial = 1;
% loop through all data
for i = 1:length(expts)
    % all files and folders in experiment folder
    trialFolders = dir([exptPath filesep expts{i} filesep 'trial*']);
    
    % loop through all trial folders in this experiment folder
    for j = 1:length(trialFolders)
        cd([trialFolders(j).folder filesep trialFolders(j).name]);
        
        trialFiles = dir('*.txt');
        
        if (length(trialFiles) == 3)
            % read in cnc data
            cncFile = fopen('cnc.txt');
            cncDat = textscan(cncFile, '%f %f %f\n', ...
                'Delimiter', ',', 'HeaderLines', 1);
            fclose(cncFile);

            trialDat(trial).cnc.t = cncDat{1};
            trialDat(trial).cnc.x = cncDat{2};
            trialDat(trial).cnc.y = cncDat{3};

            % read in cam data
            camFile = fopen('cam.txt');
            camDat = textscan(camFile, '%f %s %f %f %f %f %f\n', ...
                'Delimiter', ',', 'HeaderLines', 1);
            fclose(camFile);

            trialDat(trial).cam.t = camDat{1};
            trialDat(trial).cam.flyPresent = strcmpi(camDat{2}, 'True');
            trialDat(trial).cam.x = camDat{3};
            trialDat(trial).cam.y = camDat{4};
            trialDat(trial).cam.ang = camDat{7};

            % read in stimulus data
            stimFile = fopen('stimuli.txt');
            stimTemp = textscan(stimFile, '%s', 'Delimiter', {':',','});
            for k = 3:2:(length(stimTemp{1})-2)
                % extract variable name
                varName = stimTemp{1}{k}(2:(end-1));
                % extract variable value
                if (contains(stimTemp{1}{k+1}, '"'))
                    varVal = stimTemp{1}{k+1}(2:(end-1));
                else
                    varVal = str2num(stimTemp{1}{k+1});
                end

                stimDat.(varName) = varVal;
            end
            trialDat(trial).stim = stimDat;

            fclose(stimFile);

            % increment counter
            trial = trial + 1;   
        end
    end
end

%% process data - convert to "fly" data structure
for i = 1:length(trialDat)
    
    % if the fly was present at all during the trial
    if (any(trialDat(i).cam.flyPresent))
        
        % camera data only for when fly present
        camT = trialDat(i).cam.t(trialDat(i).cam.flyPresent);
        camX = trialDat(i).cam.x(trialDat(i).cam.flyPresent);
        camY = trialDat(i).cam.y(trialDat(i).cam.flyPresent);
        camAng = trialDat(i).cam.ang(trialDat(i).cam.flyPresent);
        
        % start and end times of trajectory
        tStart = max([camT(1), trialDat(i).cnc.t(1)]);
        tEnd = min([camT(end), trialDat(i).cnc.t(end)]);
        
        % get fly's interpolated time, x position, y position, and angle
        fly(i).t = tStart:TIME_RES:tEnd;
        fly(i).tZeroed = fly(i).t - fly(i).t(1);
        
        % fly's x, y position is sum of camera and and cnc position
        fly(i).x = interp1(camT, camX, fly(i).t) + ...
            interp1(trialDat(i).cnc.t, trialDat(i).cnc.x, fly(i).t);
        fly(i).y = interp1(camT, camY, fly(i).t) + ...
            interp1(trialDat(i).cnc.t, trialDat(i).cnc.y, fly(i).t);
        fly(i).rawAngs = interp1(camT, camAng, fly(i).t);
        
        % unwrap angles        
        jumpThresh = 150;
        unWrAngs = [];
        angOffset = 0;

        for j = 1:length(fly(i).rawAngs)
            curAng = fly(i).rawAngs(j);
            if (j == 1)
                unWrAngs = [unWrAngs curAng];
                continue;
            end

            curAng = curAng + angOffset;

            if ((unWrAngs(end) - curAng) > jumpThresh)
                angOffset = angOffset + 180;
                curAng = curAng + 180;
            elseif ((curAng - unWrAngs(end)) > jumpThresh)
                angOffset = angOffset - 180;
                curAng = curAng - 180;
            end

            unWrAngs = [unWrAngs curAng];
        end
        
        % correct angles from ellipse fitting definition, camera flip
        unWrAngs = 90 - unWrAngs;
        unWrAngs = -1 * unWrAngs;
        
        fly(i).unWrAng = unWrAngs;

        
        % smoothing
        if (SIGMA ~= 0)
            fly(i).x = smoothdata(fly(i).x, 'gaussian', SIGMA);
            fly(i).y = smoothdata(fly(i).y, 'gaussian', SIGMA);
            fly(i).unWrAng = smoothdata(fly(i).unWrAng, 'gaussian', SIGMA);
        end
        
        fly(i).wrAng = wrapTo360(unWrAngs); % wrapped angles
        
        % compute velocities
        xVel = gradient(fly(i).x) .* (1/TIME_RES);
        yVel = gradient(fly(i).y) .* (1/TIME_RES);
        
        % fly's angular velocity, turns to fly's own right are positive
        fly(i).angVel = gradient(fly(i).unWrAng) .* (1/TIME_RES);
        
        % translational velocity in any direction
        fly(i).transVel = sqrt(xVel.^2 + yVel.^2);
        
        % forward velocity
        velAng = atan2d(yVel, xVel);
        angDiff = velAng - fly(i).wrAng; 
        % fly's forward velocity, forward is positive
        fly(i).fwdVel = fly(i).transVel .* cosd(angDiff);
        
        % fly's lateral velocity, movement to fly's own right is positive
        fly(i).latVel = fly(i).transVel .* sind(angDiff);  
        
        % distance covered
        xCovered = sum(abs(diff(fly(i).x)));
        yCovered = sum(abs(diff(fly(i).y)));
        fly(i).distCovered = xCovered + yCovered;
        fly(i).angCovered = sum(abs(diff(fly(i).unWrAng)));
        
        % duration of trial
        fly(i).duration = fly(i).t(end) - fly(i).t(1);
        
        % assign fly types based on stimulus
        switch trialDat(i).stim.name
            case 'ConstantBackground'
                switch trialDat(i).stim.background
                    case 0
                        fly(i).type = 'dark';
                    case 0.5
                        fly(i).type = 'gray';
                    case 1
                        fly(i).type = 'bright';
                end
            case 'SineGrating'
                switch trialDat(i).stim.angle
                    case 0
                        fly(i).type = 'vertical';
                    case 90
                        fly(i).type = 'horizontal';
                end
            case 'RandomGrid'
                fly(i).type = 'checker';
            otherwise
                fly(i).type = [];
        end 
    end
end

%% save data

% save processed data
save([exptPath filesep 'pDatAll.mat'], 'fly', '-v7.3');
% save raw data, separately
save([exptPath filesep 'rawDatAll.mat'], 'trialDat', '-v7.3');
