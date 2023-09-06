% getAlignedLegStepsContfromBouts.m
%
% Helper function for saveContLegStepParamCond_bouts() that takes in 
%  indices of yaw velocity peaks (peak, start, and end) and returns the 
%  continuous leg step parameter values aligned to the velocity peak
% To allow averaging over different bouts and the inconsistent frame rate,
%  interpolate to specified frame rate
%
% INPUTS:
%   whichParam - string for name of legStepsCont parameter to get aligned
%       value for
%   legStepsCont - struct of continuous leg step parameter values
%   validDur - time in seconds on each side of peak to consider
%   interpFramerate - frame rate to interpolate output continuous leg
%       parameters
%   peakInd - indices of frames of bout peaks
%   flipLR - boolean for whether to flip the leg assignments left/right
%   matchedLegInd - indices for pairing of left and right legs
%   inv - boolean for whether to sign invert the parameter value
%
% OUTPUTS:
%   outContVal - numTimePts x numLegs x numBouts matrix for value of 
%       specified legStepCont parameter, for bouts, aligned to yaw velocity 
%       peak
%   t - time vector for output vals, interpolated
%
% CREATED: 9/5/23 - HHY
%
% UPDATED:
%   9/5/23 - HHY
%
function [outContVal, t] = getAlignedLegStepsContfromBouts(...
    whichParam, legStepsCont, validDur, interpFramerate, ...
    peakInd, flipLR, matchedLegInd, inv)

    % get number time pts to either side of peak, after interpolation
    numFrames = ceil(validDur * interpFramerate);

    % number of bouts
    numBouts = length(peakInd);
    % number of tracked pts
    numLegs = length(legStepsCont.legIDs.ind);
    % interframe interval, for interpolation
    ifi = 1/interpFramerate;

    % preallocate
    outContVal = nan(numFrames * 2 + 1, numLegs, numBouts);


    % loop through all bouts
    for i = 1:numBouts
        % zero time vector, so that yaw velocity peaks are at t = 0
        tOrig = legStepsCont.t - legStepsCont.t(peakInd(i));

        % get time vector for interpolation - only consider +/- validDur
        %  around peak
        % doing it this way keeps 0 at 0
        newTDur = (numFrames / interpFramerate);


        newTHalf1 = 0:ifi:newTDur;
        newTHalf1 = fliplr(newTHalf1) * -1;

        newTHalf2 = 0:ifi:newTDur;

        newT = [newTHalf1 newTHalf2(2:end)];

        % get interpolated continuous step parameter values
        interpStepVal = interp1(tOrig, legStepsCont.(whichParam), ...
            newT, 'spline');

        interpStepValFlip = interpStepVal;

        % if flipping legs L/R
        if (flipLR)
            interpStepValFlip(:,matchedLegInd(:,1)) = ...
                interpStepVal(:,matchedLegInd(:,2));
            interpStepValFlip(:,matchedLegInd(:,2)) = ...
                interpStepVal(:,matchedLegInd(:,1));
        end

        % if inverting the parameter value
        if (inv)
            interpStepValFlip = interpStepValFlip * -1;
        end

        % add to output matrix
        outContVal(:,:,i) = interpStepValFlip;
    end

    % get time vector for output matrix
    tHalf1 = fliplr((0:numFrames) * ifi * -1);
    tHalf2 = (1:numFrames) * ifi;

    t = [tHalf1 tHalf2]';
end