% preprocessFreeWalkLegTrack.m
%
% Function that takes in full path to .trk file and returns preprocessed
%  leg tracking data (leg position and velocity) in struct.
% .trk file is output from APT
%
% Modification of preprocessLegTrack on tethered flies
%
% INPUTS:
%   trkFilepath - full path to .trk file 
%   frameTimes - time at which each leg video frame occurred
%   startFrame - start frame for time period of video to consider
%   endFrame - end frame for time period of video to consider
%   cncX - cnc x positions, in meters
%   cncY - cnc y positions, in meters
%   refPts - struct of indicies corresponding to reference points
%       headPtInd - index of head of fly
%       thxRPtInd - index of thorax right point
%       thxLPtInd - index of thorax left point
%       thxMidInd - index of midpoint of back of thorax
%       abdPtInd - index of abdomen tip of fly
%   smoParams - struct of parameters for smoothing leg velocity
%   imgParams - struct of video params
%       PX_PER_METER - pixels per meter conversion
%       xPx - number of pixels in x dimension (width)
%       yPx - number of pixels in y dimension (height)
%
% OUTPUTS:
%   legTrack - struct with processed data, fields:
%       trkFilepath - full path to .trk file
%       trkFilename - name of .trk file
%       legX - raw leg positions in X (front-back axis), straight from .trk
%           but only between start and end frames
%       legY - raw leg positions in Y (front-back axis), straight from .trk
%           but only between start and end frames
%       numFrames - number of leg tracking frames
%       numPoints - number of points being tracked
%       t - time vector for leg position/velocity data, time of each leg
%           video frame
%       srnLegX - leg positions in X after alignment and normalization
%       srnLegY - leg positions in Y after alignment and normalization
%       aLegX - absolute leg positions in X, in meters
%       aLegY - absolute leg positions in Y, in meters
%       srnfLegX - leg positions in X after alignment and normalization and
%           with outliers removed by median filtering
%       srnfLegY - leg positions in Y after alignment and normalization
%           with outliers removed by median filtering
%       legXVel - instantaneous leg velocity in X, with Gaussian process
%           smoothing
%       legYVel - instantaneous leg velocity in Y, with Gaussian process
%           smoothing
%       refPts - struct of parameters defining reference tracked points
%       smoParams - struct of parameters for smoothing leg velocity
%
% CREATED: 4/8/23 - HHY
%
% UPDATED:
%   4/11/23 - HHY
%
function legTrack = preprocessFreeWalkLegTrack(trkFilepath, frameTimes, ...
    startFrame, endFrame, cncX, cncY, refPts, smoParams, imgParams)

    % load .trk file; get leg positions
    [legX, legY] = loadTrkFile(trkFilepath);

    % tracking only b/w start and end frames
    legTrack.legX = legX(startFrame:endFrame,:);
    legTrack.legY = legY(startFrame:endFrame,:);

    % some parameters about .trk file
    % number of frames in this trial
    legTrack.numFrames = size(legTrack.legX, 1); 
    % number of tracked points
    legTrack.numPoints = size(legTrack.legX, 2);
    
    % frame times
    legTrack.t = frameTimes;

    % get leg positions after alignment to fly midpoint and normalization 
    %  to body length units
    [legTrack.srnLegX, legTrack.srnLegY, legTrack.aLegX, ...
        legTrack.aLegY, legTrack.flyRawAng] = ...
        getLegBodyPosFromTrk(legTrack.legX, legTrack.legY, cncX, cncY,...
        refPts, imgParams);

    % remove outliers from leg position, ones that deviate too far from
    %  median value
    for i = 1:legTrack.numPoints
        legTrack.srnfLegX(:,i) = medFiltRmvOutliers(...
            legTrack.srnLegX(:,i), smoParams.medFiltDeg, ...
            smoParams.percntDev, smoParams.maxMinWin);
        legTrack.srnfLegY(:,i) = medFiltRmvOutliers(...
            legTrack.srnLegY(:,i), smoParams.medFiltDeg, ...
            smoParams.percntDev, smoParams.maxMinWin);
    end

    % get leg velocities, with Gaussian process smoothing on position
    legTrack.legXVel = findLegVel(legTrack.srnfLegX, smoParams);
    legTrack.legYVel = findLegVel(legTrack.srnfLegY, smoParams);
    
    % copy over parameter structs
    legTrack.refPts = refPts;
    legTrack.smoParams = smoParams;

end