% getAlignedLegXYfromBouts.m
%
% Helper function for saveLegStepParamCond_bouts() that takes in indices of
%  yaw velocity peaks (peak, start, and end) and returns the leg X and Y
%  positions aligned to the velocity peak
% To allow averaging over different bouts and the inconsistent frame rate,
%  interpolate to specified frame rate
%
% INPUTS:
%   legTrack - struct of leg tracking info, including leg X and Y positions
%   legXYParams - struct of parameters, directly from
%     saveLegStepParamCond_bouts()
%       maxDuration - time in seconds to consider on each side 
%       interpFrameRate - frame rate to interpolate to, in Hz
%   peakInd - indices of frames of bout peaks
%   boutStartInd - indices of bout starts
%   boutEndInd - indices of bout ends
%
% OUTPUTS:
%   boutLegX - numTimePts x numTrackedPts x numBouts matrix for leg X
%       positions, for bouts, aligned to yaw velocity peak
%   boutLegY - like boutLegX, but for leg Y positions
%   t - time vector for positions, interpolated
%
% CREATED: 5/12/23 - HHY
%
% UPDATED:
%   5/12/23 - HHY
%
function [boutLegX, boutLegY, t] = getAlignedLegXYfromBouts(legTrack, ...
    legXYParams, peakInd, boutStartInd, boutEndInd)

    % get number time pts to either side of peak, after interpolation
    maxNumFrames = floor(legXYParams.maxDuration * ...
        legXYParams.interpFrameRate);
    % number of bouts
    numBouts = length(peakInd);
    % number of tracked pts
    numTrkPts = legTrack.numPoints;

    % preallocate
    boutLegX = nan(maxNumFrames * 2 + 1, numTrkPts, numBouts);
    boutLegY = nan(maxNumFrames * 2 + 1, numTrkPts, numBouts);


    % loop through all bouts
    for i = 1:numBouts
        % zero time vector, so that yaw velocity peaks are at t = 0
        tOrig = legTrack.t - legTrack.t(peakInd(i));

        % get time for bout start and end
        boutStartT = legTrack.t(boutStartInd(i))- legTrack.t(peakInd(i));
        boutEndT = legTrack.t(boutEndInd(i))- legTrack.t(peakInd(i));

        % get time vector for interpolation - only consider +/- maxDuration
        %  around peak
        % doing it this way keeps 0 at 0
        newTDur = (maxNumFrames / legXYParams.interpFrameRate);

        ifi = 1/legXYParams.interpFrameRate;

        newTHalf1 = 0:ifi:newTDur;
        newTHalf1 = fliplr(newTHalf1) * -1;

        newTHalf2 = 0:ifi:newTDur;

        newT = [newTHalf1 newTHalf2(2:end)];

        % get interpolated leg X and Y position
        interpX = interp1(tOrig, legTrack.srnfLegX, newT,'spline');
        interpY = interp1(tOrig, legTrack.srnfLegY, newT,'spline');

        % filter for bout start and end
        boutStartLog = newT < boutStartT;
        interpX(boutStartLog,:) = nan;
        interpY(boutStartLog,:) = nan;

        boutEndLog = newT > boutEndT;
        interpX(boutEndLog,:) = nan;
        interpY(boutEndLog,:) = nan;

        % add to output matrix
        boutLegX(:,:,i) = interpX;
        boutLegY(:,:,i) = interpY;

    end
    
end