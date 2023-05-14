% computeFreeWalkStepParameters.m
%
% Function that takes step indices [start, mid, end] and computes step
%  parameters: step lengths, step directions, step durations, step speeds,
%  step X velocities, step Y velocities
%
% Adapated from computeStepParameters(), but only computes step parameters
%  relevant to free walking data
%
% INPUTS:
%   legSteps - struct of leg step data, from getLegReversals() and 
%     minMaxPosToSteps() with fields:
%       maxIndsAll
%       minIndsAll
%       maxWhichLeg
%       minWhichLeg
%       stepInds
%       stepWhichLeg
%   legTrack - struct of leg tracking data, output of 
%     loadFreeWalkTrk2PData()
%   bodytraj - struct of body trajectory data, output of
%     preprocessBodyTraj()
%   
% OUTPUTS:
%   legSteps - struct of leg step data, updated with additional fields.
%     Each new field is n x 2 matrix of steps x (val start to mid), (val
%     mid to end), except stepT, which is n x 3 (t at start, mid, end)
%       stepLengths - step length in xy plane
%       stepXLengths - step length in x direction
%       stepYLengths - step length in y direction
%       stepDirections - step direction, degrees, in xy plane
%       stepDurations - step duration (in sec)
%       stepSpeeds - step speed in xy plane (body lengths/sec)
%       stepVelX - step velocity in x direction (body lengths/sec)
%       stepVelY - step velocity in y direction (body lengths/sec)
%       stepAEPX - step anterior extreme position in X
%       stepAEPY - step Y position for AEP as defined in X
%       stepPEPX - step posterior extreme position in X
%       stepPEPY - step Y position for PEP as defined in X
%       stepAbsDisp - step absolute displacement in arena frame (m)
%       stepAbsTransSpd - step absolute translational speed in arena frame
%           (m/sec)
%       stepT - time at start, mid, and end points of step (sec)
%       stepBodyFwd - body forward velocity during step (m/sec)
%       stepBodyLat - body lateral velocity during step (m/sec)
%       stepBodyYaw - body lateral velocity during step (deg/sec)
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%
function legSteps = computeFreeWalkStepParameters(legSteps, legTrack, ...
    bodytraj)

    % preallocate
    stepLengths = zeros(size(legSteps.stepInds,1), 2); 
    stepXLengths = zeros(size(legSteps.stepInds,1), 2);  
    stepYLengths = zeros(size(legSteps.stepInds,1), 2);
    stepDirections = zeros(size(legSteps.stepInds,1), 2);
    stepDurations = zeros(size(legSteps.stepInds,1), 2);
    stepSpeeds = zeros(size(legSteps.stepInds,1), 2);
    stepVelX = zeros(size(legSteps.stepInds,1), 2);
    stepVelY = zeros(size(legSteps.stepInds,1), 2);
    stepAEPX = zeros(size(legSteps.stepInds,1),2);
    stepAEPY = zeros(size(legSteps.stepInds,1),2);
    stepPEPX = zeros(size(legSteps.stepInds,1),2);
    stepPEPY = zeros(size(legSteps.stepInds,1),2);
    stepAbsDisp = zeros(size(legSteps.stepInds,1),2);
    stepAbsTransSpd = zeros(size(legSteps.stepInds,1),2);
    stepT = zeros(size(legSteps.stepInds,1),3);
    stepBodyFwd = zeros(size(legSteps.stepInds,1),2);
    stepBodyYaw = zeros(size(legSteps.stepInds,1),2);
    stepBodyLat = zeros(size(legSteps.stepInds,1),2);
    
    % loop through all steps, compute parameters for each half step
    for i = 1:size(legSteps.stepInds, 1)
        % get X position for step start
        stepStartX = legTrack.srnfLegX(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step start
        stepStartY = legTrack.srnfLegY(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get X position for step mid
        stepMidX = legTrack.srnfLegX(legSteps.stepInds(i,2), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step mid
        stepMidY = legTrack.srnfLegY(legSteps.stepInds(i,2), ...
            legSteps.stepWhichLeg(i));
        % get X position for step end
        stepEndX = legTrack.srnfLegX(legSteps.stepInds(i,3), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step end
        stepEndY = legTrack.srnfLegY(legSteps.stepInds(i,3), ...
            legSteps.stepWhichLeg(i));
        
        % get absolute X position for step start
        stepStartAbsX = legTrack.aLegX(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get absolute Y position for step start
        stepStartAbsY = legTrack.aLegY(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get absolute X position for step mid
        stepMidAbsX = legTrack.aLegX(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get absolute Y position for step mid
        stepMidAbsY = legTrack.aLegY(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get absolute X position for step end
        stepEndAbsX = legTrack.aLegX(legSteps.stepInds(i,1),...
            legSteps.stepWhichLeg(i));
        % get absolute Y position for step end
        stepEndAbsY = legTrack.aLegY(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        

        % get length of first half step (start to mid)
        stepLengths(i,1) = distBtw2Pts(stepStartX, stepStartY, ...
            stepMidX, stepMidY);
        % get length of second half step (mid to end)
        stepLengths(i,2) = distBtw2Pts(stepMidX, stepMidY, stepEndX, ...
            stepEndY);
        
        % step X lengths
        stepXLengths(i,1) = stepMidX - stepStartX;
        stepXLengths(i,2) = stepEndX - stepMidX;

        % step Y lengths
        stepYLengths(i,1) = stepMidY - stepStartY;
        stepYLengths(i,2) = stepEndY - stepMidY;

        % get absolute displacement for first half step (start to mid)
        stepAbsDisp(i,1) = distBtw2Pts(stepStartAbsX, stepStartAbsY, ...
            stepMidAbsX, stepMidAbsY);
        % second half step
        stepAbsDisp(i,2) = distBtw2Pts(stepMidAbsX, stepMidAbsY, ...
            stepEndAbsX, stepEndAbsY);

        % step AEP
        % get whether start or mid position is most anterior (min position)
        [~, aepInd] = min([stepStartX, stepMidX]);
        % match AEP X to AEP Y; PEP is one that's not AEP
        if (aepInd == 1)
            stepAEPX(i,1) = stepStartX;
            stepAEPY(i,1) = stepStartY;
            stepPEPX(i,1) = stepMidX;
            stepPEPY(i,1) = stepMidY;
        else
            stepAEPX(i,1) = stepMidX;
            stepAEPY(i,1) = stepMidY;
            stepPEPX(i,1) = stepStartX;
            stepPEPY(i,1) = stepStartY;
        end
        % same for second half step
        [~, aepInd] = min([stepMidX, stepEndX]);
        if (aepInd == 1)
            stepAEPX(i,2) = stepMidX;
            stepAEPY(i,2) = stepMidY;
            stepPEPX(i,2) = stepEndX;
            stepPEPY(i,2) = stepEndY;
        else
            stepAEPX(i,2) = stepEndY;
            stepAEPY(i,2) = stepEndY;
            stepPEPX(i,2) = stepMidX;
            stepPEPY(i,2) = stepEndY;
        end

        % get direction of first half step (start to mid)
        stepDirections(i,1) = findAngle2Pts(stepStartX, stepStartY, ...
            stepMidX, stepMidY);
        % get direction of second half step (mid to end)
        stepDirections(i,2) = findAngle2Pts(stepMidX, stepMidY, ...
            stepEndX, stepEndY);

        % get time for step start
        stepStartT = legTrack.t(legSteps.stepInds(i,1));
        % get time for step mid
        stepMidT = legTrack.t(legSteps.stepInds(i,2));
        % get time for step end
        stepEndT = legTrack.t(legSteps.stepInds(i,3));

        % step time
        stepT(i,:) = [stepStartT, stepMidT, stepEndT];

        % get duration of first half step (start to mid)
        stepDurations(i,1) = stepMidT - stepStartT;
        % get duration of second half step (mid to end)
        stepDurations(i,2) = stepEndT - stepMidT;

        % get velocity in xy plane of first half step (start to mid)
        stepSpeeds(i,1) = stepLengths(i,1) / stepDurations(i,1);
        % get velocity in xy plane of second half step (mid to end)
        stepSpeeds(i,2) = stepLengths(i,2) / stepDurations(i,2);

        % get velocity in x for first half of step (start to mid)
        stepVelX(i,1) = (stepMidX - stepStartX) / stepDurations(i,1);
        % get velocity in x for second half step (mid to end)
        stepVelX(i,2) = (stepEndX - stepMidX) / stepDurations(i,2);

        % get velocity in y for first half of step (start to mid)
        stepVelY(i,1) = (stepMidY - stepStartY) / stepDurations(i,1);
        % get velocity in y for second half step (mid to end)
        stepVelY(i,2) = (stepEndY - stepMidY) / stepDurations(i,2); 

        % get translational speed in arena reference frame
        % first half step
        stepAbsTransSpd(i,1) = stepAbsDisp(i,1) / stepDurations(i,1);
        % second half step
        stepAbsTransSpd(i,2) = stepAbsDisp(i,2) / stepDurations(i,2);
        
        % get net displacement of forward position during first half step
        stepBodyFwdDispl1 = trapz(...
            bodytraj.tZeroed(legSteps.stepInds(i,1):legSteps.stepInds(i,2)), ...
            bodytraj.fwdVelSmoM((legSteps.stepInds(i,1):legSteps.stepInds(i,2))));
        % net displacement of forward position during second half step
        stepBodyFwdDispl2 = trapz(...
            bodytraj.tZeroed(legSteps.stepInds(i,2):legSteps.stepInds(i,3)), ...
            bodytraj.fwdVelSmoM((legSteps.stepInds(i,2):legSteps.stepInds(i,3))));

        % fwd velocity, first half step (divide displacement by duration)
        stepBodyFwd(i,1) = stepBodyFwdDispl1 / stepDurations(i,1);
        % second half step
        stepBodyFwd(i,2) = stepBodyFwdDispl2 / stepDurations(i,2);

        % angular position for step start point
        stepStartYawPos = bodytraj.unWrAngSmoM(legSteps.stepInds(i,1));
        % FicTrac yaw position for step mid point
        stepMidYawPos = bodytraj.unWrAngSmoM(legSteps.stepInds(i,2));
        % FicTrac yaw position for step end point
        stepEndYawPos = bodytraj.unWrAngSmoM(legSteps.stepInds(i,3));

        % FicTrac yaw velocity, first half step
        stepBodyYaw(i,1) = (stepMidYawPos-stepStartYawPos) / ...
            stepDurations(i,1);
        % second half step
        stepBodyYaw(i,2) = (stepEndYawPos - stepMidYawPos) / ...
            stepDurations(i,2);

        % get net displacement of lateral position during first half step
        stepBodyLatDispl1 = trapz(...
            bodytraj.tZeroed(legSteps.stepInds(i,1):legSteps.stepInds(i,2)), ...
            bodytraj.latVelSmoM((legSteps.stepInds(i,1):legSteps.stepInds(i,2))));
        % net displacement of lateral position during second half step
        stepBodyLatDispl2 = trapz(...
            bodytraj.tZeroed(legSteps.stepInds(i,2):legSteps.stepInds(i,3)), ...
            bodytraj.latVelSmoM((legSteps.stepInds(i,2):legSteps.stepInds(i,3))));

        % lat velocity, first half step (divide displacement by duration)
        stepBodyLat(i,1) = stepBodyLatDispl1 / stepDurations(i,1);
        % second half step
        stepBodyLat(i,2) = stepBodyLatDispl2 / stepDurations(i,2);

    end
    
    % add these parameters to legSteps struct
    legSteps.stepLengths = stepLengths;
    legSteps.stepXLengths = stepXLengths;
    legSteps.stepYLengths = stepYLengths;
    legSteps.stepDirections = stepDirections;
    legSteps.stepDurations = stepDurations;
    legSteps.stepSpeeds = stepSpeeds;
    legSteps.stepVelX = stepVelX;
    legSteps.stepVelY = stepVelY;
    legSteps.stepAEPX = stepAEPX;
    legSteps.stepAEPY = stepAEPY;
    legSteps.stepPEPX = stepPEPX;
    legSteps.stepPEPY = stepPEPY;
    legSteps.stepAbsDisp = stepAbsDisp;
    legSteps.stepAbsTransSpd = stepAbsTransSpd;
    legSteps.stepT = stepT;
    legSteps.stepBodyFwd = stepBodyFwd;
    legSteps.stepBodyLat = stepBodyLat;
    legSteps.stepBodyYaw = stepBodyYaw;
end