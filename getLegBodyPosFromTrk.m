% getLegBodyPosFromTrk.m
%
% Function that takes in raw leg X and Y coordinates (output from
%  loadTrkFile) and returns coordinates that are shifted and rotated to
%  align all the frames to the center point of the fly (point at middle of
%  back of thorax) and normalized so the units are body lengths (instead of
%  pixels). Also returns coordinates in arena location units (meters).
%  Returns fly's heading angle. 
%
% Adapted from shiftRotateNormalizeLegPos()
%
% NOTE: currently doesn't account for absolute heading angle frame of 2022
%   data (b/c of bug in acquisition timing, don't feel like resolving right
%   now - 4/11/23). Should be flyRawAng - 90 to align.
%
% INPUTS:
%   legX - X coordinates of leg positions (frames x points), output from 
%       loadTrkFile
%   legY - Y coordinates of leg positions (frames x points), output from 
%       loadTrkFile
%   cncX - cnc x positions, in meters
%   cncY - cnc y positions, in meters
%   startFrame - start frame for time period of video to consider
%   endFrame - end frame for time period of video to consider
%   refPts - struct of indicies corresponding to reference points
%       headPtInd - index of head of fly
%       thxRPtInd - index of thorax right point
%       thxLPtInd - index of thorax left point
%       thxMidInd - index of midpoint of back of thorax
%       abdPtInd - index of abdomen tip of fly
%   imgParams - struct of video params
%       PX_PER_METER - pixels per meter conversion
%       xPx - number of pixels in x dimension (width)
%       yPx - number of pixels in y dimension (height)
%
%
% OUTPUTS:
%   srnLegX - X coordinates of leg positions (frames x points), after
%       shifting and rotating to align and normalizing to body lengths 
%       units
%   srnLegY - Y coordinates of leg positions (frames x points), after
%       shifting and rotating to align and normalizing to body lengths 
%       units
%   aLegX - X coordinates of leg positions (frames x points), in absolute 
%       arena location (meters)
%   aLegY - Y coordinates of leg positions (frames x points), in absolute 
%       arena location (meters)
%   flyRawAng - fly's heading (at each frame), computed as angle of line
%       fit to 
%
% CREATED: 4/9/23 - HHY
%
% UPDATED:
%   4/9/23 - HHY
%   4/11/23 - HHY - checked x,y coordinates, angle; aligns b/w body and leg
%       tracking (for 2018 data, angle can be 180 deg offset b/c
%       rostral-caudal ambiguity in body tracking)
%
function [srnLegX, srnLegY, aLegX, aLegY, flyRawAng] = ...
    getLegBodyPosFromTrk(legX, legY, cncX, cncY, refPts, imgParams)

    % get positions of tracked points, relative to fly coordinate frame

    % get midpoint between thorax left and right points, front of thorax
    thxFrntX = (legX(:,refPts.thxRPtInd) + legX(:,refPts.thxLPtInd))/2;
    thxFrntY= (legY(:,refPts.thxRPtInd) + legY(:,refPts.thxLPtInd))/2;
    
    % fit line to head, thx front, thx back, abd

    % preallocate vector for saving coefficients of line fit
    bodyLineCoeffs = zeros(size(legX,1), 2);

    % indices of points to fit line
    lineFitInd = [refPts.headPtInd, refPts.thxMidInd, refPts.abdPtInd];

    % get coefficients of line for each video frame; first is slope, second
    %  is y-intercept
    for i = 1:size(legX, 1)
        xPts = legX(i, lineFitInd);
        yPts = legY(i, lineFitInd);

        % add in thx front point
        xPts = [xPts thxFrntX(i)];
        yPts = [yPts thxFrntY(i)];

        bodyLineCoeffs(i,:) = polyfit(xPts, yPts, 1);
    end

    % get projection of back of thorax point onto body line (define 
    %  this as midpoint of fly)
    [midXProj, midYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
        bodyLineCoeffs(:,2), legX(:,refPts.thxMidInd), ...
        legY(:,refPts.thxMidInd));

    % project head point and abdomen point onto body line, get distance
    %  between them to get body length
    [headXProj, headYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
        bodyLineCoeffs(:,2), legX(:,refPts.headPtInd), ...
        legY(:,refPts.headPtInd));
    [abdXProj, abdYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
        bodyLineCoeffs(:,2), legX(:,refPts.abdPtInd), ...
        legY(:,refPts.abdPtInd));

    bodyLen = distBtw2Pts(headXProj, headYProj, abdXProj, abdYProj);

    % shift all points so that they're relative to the midpoint, which is
    %  now defined as (0,0)
    % preallocate
    legXShift = zeros(size(legX));
    legYShift = zeros(size(legY));
    % shift each point
    for i = 1:size(legX,2)
        legXShift(:,i) = legX(:,i) - midXProj;
        legYShift(:,i) = legY(:,i) - midYProj;
    end
    thxFrntXShift = thxFrntX - midXProj;
    thxFrntYShift = thxFrntY - midYProj;

    % get body line coefficients, after shift
    bodyLineShiftCoeffs = zeros(size(legXShift,1),2);

    for i = 1:size(legXShift, 1)
        xPts = legXShift(i, lineFitInd);
        yPts = legYShift(i, lineFitInd);

        % add in thx front point
        xPts = [xPts thxFrntXShift(i)];
        yPts = [yPts thxFrntYShift(i)];

        bodyLineShiftCoeffs(i,:) = polyfit(xPts, yPts, 1);
    end

    % for whole time series, define the body length and the body angle
    % for each frame, rotate around the midpoint by the body angle so x and
    %  y axes align with forward and lateral axes of fly
    % normalize units to body length (from pixels)

    % body length as median of body lengths throughout trial
    flyBodyLength = median(bodyLen);

    % tangent of body angle = slope of body line, 
    %  post shift (y-int now ~0)
    flyBodyAng = atand(bodyLineShiftCoeffs(:,1));

    % rotate all points about midpt, now (0,0), by inverse of fly body 
    %  angle
    % preallocate
    legXShiftRot = zeros(size(legXShift));
    legYShiftRot = zeros(size(legYShift));
    % loop through each marked point and rotate
    for i = 1:size(legXShift,2)
        % each frame
        for j = 1:size(legXShift,1)
            [legXShiftRot(j,i), legYShiftRot(j,i)] = ...
                rotatePts2D(legXShift(j,i), legYShift(j,i), -1*flyBodyAng(j));
        end
    end

    % rotation doesn't account for which side is head and which is abd
    % align such that head x < abd x
    % loop through each labeled point
    legXShiftRot2 = legXShiftRot;
    legYShiftRot2 = legYShiftRot;
    for i = 1:size(legXShiftRot,2)
        % loop through each frame
        for j = 1:size(legXShiftRot,1)
            % check if head x is > abd x (need to rotate this frame if yes)
            if (legXShiftRot(j,refPts.headPtInd) > ...
                    legXShiftRot(j,refPts.abdPtInd))
                % rotate 180 deg
                [legXShiftRot2(j,i), legYShiftRot2(j,i)] = ...
                    rotatePts2D(legXShiftRot(j,i), legYShiftRot(j,i), 180);
                % also, add 180 deg to flyBodyAng
                % to align to 2018 ref frame (mostly, that heading can be
                %  180 deg off b/c unknown rostral-caudal directionality)
                flyBodyAng(j) = flyBodyAng(j) + 180;
            end
        end
    end

    % wrap flyBodyAng to 360 deg to get flyRawAng
    flyRawAng = wrapTo360(flyBodyAng);

    % normalize to body length
    srnLegX = legXShiftRot2 / flyBodyLength;
    srnLegY = legYShiftRot2 / flyBodyLength;

    % multiply srnLegY by -1 so that fly's right side has positive
    %  coordinates and left side has negative coordinates
    srnLegY = -1 * srnLegY;


    % get absolute positions of tracked points, in meters
    aLegX = cncX - ((legX - (imgParams.xPx/2)) / imgParams.PX_PER_METER);
    aLegY = cncY - ((legY - (imgParams.yPx/2)) / imgParams.PX_PER_METER);

end