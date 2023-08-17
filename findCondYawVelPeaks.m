% findCondYawVelPeaks.m
%
% Helper function for saveLegStepParamCond_bouts() that takes in bodytraj,
%  cond, and moveNotMove and returns indices (of bodytraj) for valid yaw
%  velocity peaks as well as start/end indices for corresponding bouts.
% Boolean for whether to extract right or left turns
%
% INPUTS:
%   bodytraj - pData output struct
%   cond - pass from saveLegStepParamCond_bouts()
%   moveNotMove - pData output struct
%   rightTurn - boolean for whether to extract right turns (false = left)
%   fwdVelCond - pass from saveLegStepParamCond_bouts()
%
% OUTPUTS:
%   yawVelPeakInd - vector of indices of bodytraj corresponding to yaw
%       velocity peaks (that meet cond criteria), length n
%   boutStartInd - vector of indices corresponding to start of bouts,
%       length n
%   boutEndInd - vector of indices corresponding to end of bouts, length n
%
% CREATED: 5/10/23 - HHY
%
% UPDATED:
%   5/10/23 - HHY
%   6/8/23 - HHY - slight bug fix to only invert velocities for left turns
%       (not speeds)
%   6/9/23 - HHY - update to allow conditioning on bout duration
%   8/16/23 - HHY - update to allow conditioning on initial and change in
%       forward velocity
%
function [yawVelPeakInd, boutStartInd, boutEndInd] = findCondYawVelPeaks(...
    bodytraj, cond, fwdVelCond, moveNotMove, rightTurn)

    % trial end exclusion
    % exclude any peaks that are within x frames of trial end
    TRIAL_END_EXCLUDE = 20;

    % check if we're extracting right or left turns
    if (rightTurn)
        angVelSmoS = bodytraj.angVelSmoS;

    else
        angVelSmoS = -1 * bodytraj.angVelSmoS;
    end
        
    % find all yaw velocity peaks, no conditioning yet
    % always operates on bodytraj.angVelSmoS
    [~, pkInds] = findpeaks(angVelSmoS);

    % remove peaks during not moving times
    pkInds = setdiff(pkInds, moveNotMove.notMoveInd,'stable');

    % remove any peaks too close to end of trial
    lastValidFrame = length(bodytraj.tZeroed) - TRIAL_END_EXCLUDE;

    pkInds(pkInds>lastValidFrame) = [];

    % convert pkInds to logical (for bodytraj indices)
    pkLog = false(size(bodytraj.tZeroed));
    pkLog(pkInds) = true;

    % find peaks that meet conditions of cond
    for i = 1:length(cond.whichParam)
        % if condition is on yaw or on lateral, invert if left turns
        % don't touch if right turn
        if (rightTurn)
            thisCond = bodytraj.(cond.whichParam{i});
        else % left turn
            if (contains(cond.whichParam{i}, 'angVel', 'IgnoreCase',true) || ...
                    contains(cond.whichParam{i}, 'latVel', 'IgnoreCase',true))
                thisCond = -1 * bodytraj.(cond.whichParam{i});
            else
                thisCond = bodytraj.(cond.whichParam{i});
            end
        end
         
        % all time points that meet condition
        thisCondLog = eval(['thisCond' cond.cond{i}]);

        % intersect all points that meet condition with peaks
        pkLog = pkLog & thisCondLog;
    end

    % convert logical back to indices
    pkInds = find(pkLog);

    % all indices where yaw velocity is less than min
    yawMinInd = find(angVelSmoS < cond.minYawThresh);

    % preallocate start and end index vectors
    pkStartInd = zeros(size(pkInds));
    pkEndInd = zeros(size(pkInds));

    % find start and end for each peak
    % if start and end can't be found (too close to edge, remove peaks)
    rmvInd = [];
    for i = 1:length(pkInds)
        thisStart = find(yawMinInd < pkInds(i), 1, 'last');
        % if no index found, set peak start to beginning of trial
        if isempty(thisStart)
            rmvInd = [rmvInd i];
        else
            % +1 because presumably next index is greater than min
            pkStartInd(i) = yawMinInd(thisStart) + 1;
        end

        thisEnd = find(yawMinInd > pkInds(i), 1, 'first');
        % if no index found, set peak end to end of trial
        if isempty(thisEnd)
            rmvInd = [rmvInd i];
        else
            % -1 because this is first that is less than min
            pkEndInd(i) = yawMinInd(thisEnd) - 1;
        end
    end
    % remove peaks too close to edge
    pkInds(rmvInd) = [];
    pkStartInd(rmvInd) = [];
    pkEndInd(rmvInd) = [];

    % check if any of the peaks share edges - look only at peak starts
    %  (will be the same as peak ends)

    % get indices into bodytraj, for not unique starts 
    [~, uniInd]  = unique(pkStartInd, 'stable');
    nonUniInd = setdiff(1:numel(pkStartInd), uniInd);
    startNonUniInd = pkStartInd(nonUniInd);

    % indices of peaks themselves, for non-unique peaks
    peakNonUniInd = pkInds(nonUniInd); 


    % if peaks share edges, keep only peak with greatest yaw velocity
    if (~isempty(startNonUniInd))
        notKeepInd = []; % initialize tracker of peaks to discard
        for i = 1:length(startNonUniInd)
            thisPeakStart = startNonUniInd(i);

            % find all peaks that share this start
            sharedInd = find(pkStartInd == thisPeakStart);

            % find peak with greatest yaw velocity
            bInd = pkInds(sharedInd);
            yawVals = angVelSmoS(bInd);

            [~, keepPkInd] = max(yawVals);
            keepPkBInd = bInd(keepPkInd);
            notKeepPkBInd = setdiff(bInd, keepPkBInd, 'stable');

            if (iscolumn(notKeepPkBInd))
                notKeepPkBInd = notKeepPkBInd';
            end

            % track all peaks to discard
            notKeepInd = [notKeepInd notKeepPkBInd];
        end

        % discard not keep peaks
        % discard for peaks, get indices of ones to keep
        [pkInds, keepInd] = setdiff(pkInds,notKeepInd,'stable');
        pkStartInd = pkStartInd(keepInd);
        pkEndInd = pkEndInd(keepInd);
    end

    % loop through all peaks and check that the bout meets the duration
    %  requirements
    % initialize to keep track of peaks to remove
    rmInd = []; 
    for i = 1:length(pkInds)

        thisBoutStartT = bodytraj.tZeroed(pkStartInd(i));
        thisBoutEndT = bodytraj.tZeroed(pkEndInd(i));

        thisBoutDur = thisBoutEndT - thisBoutStartT;

        % if this bout duration is too short or too long, flag this index
        %  for deletion
        if((thisBoutDur < cond.turnDur(1)) || ...
                (thisBoutDur > cond.turnDur(2)))
            rmInd = [rmInd i];
        end
    end
    % remove any bouts that don't meet the duration criteria
    pkInds(rmInd) = [];
    pkStartInd(rmInd) = [];
    pkEndInd(rmInd) = [];

    % loop through all peaks and check that the bout meets the forward
    %  velocity requirements
    % initialize to keep track of peaks to remove
    rmInd = [];
    for i = 1:length(pkInds)
        % forward velocity at start for this bout
        thisBoutInitFwdVel = bodytraj.fwdVelSmoS(pkStartInd(i));
        % forward velocity at yaw peak for this bout
        thisBoutPeakFwdVel = bodytraj.fwdVelSmoS(pkInds(i));

        % change in forward velocity
        thisBoutChangeFwdVel = thisBoutPeakFwdVel - thisBoutInitFwdVel;

        % if this bout doesn't meet forward velocity requirements, flag
        %  this index for deletion
        if ((thisBoutInitFwdVel < fwdVelCond.initVel(1)) || ...
                (thisBoutInitFwdVel > fwdVelCond.initVel(2)) || ...
                (thisBoutChangeFwdVel < fwdVelCond.change(1)) || ...
                (thisBoutChangeFwdVel > fwdVelCond.change(2)))
            rmInd = [rmInd i];
        end
    end
    % remove any bouts that don't meet the forward velocity criteria
    pkInds(rmInd) = [];
    pkStartInd(rmInd) = [];
    pkEndInd(rmInd) = [];


    % outputs
    yawVelPeakInd = pkInds;
    boutStartInd = pkStartInd;
    boutEndInd = pkEndInd;
end