% callFreeWalkSwingStanceSteps.m
%
% Function that takes in steps (defined by start, mid, end indices),
%  step parameters (calculated in computeFreeWalkStepParameters())
% Calls swing/stance based on whether leg is moving in the reference frame
%  of the arena during that half step (swing) or the leg is stationary in 
%  the refrence frame of the arena during that half step (stance). Assign
%  half step with lower translational speed in reference frame of 
%  arena to stance (not based on thresholding)
%
% INPUTS:
%   legSteps - struct of leg step data, output from processLegTrack(), with
%       fields:
%       maxIndsAll
%       minIndsAll
%       maxWhichLeg
%       minWhichLeg
%       userSelVal
%       stepInds
%       stepWhichLeg
%       stepLengths
%       stepXLengths
%       stepYLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%       stepAEPX
%       stepAEPY
%       stepPEPX
%       stepPEPY
%       stepAbsDisp
%       stepAbsTransSpd
%       stepBodyFwd
%       stepBodyLat
%       stepBodyYaw
%
% OUTPUTS:
%   legSteps - struct of leg step data, updated with additional fields:
%       stepSwingStance - n x 2 matrix for swing/stance calls for each half
%           step; 1 for stance, -1 for swing
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%
function legSteps = callFreeWalkSwingStanceSteps(legSteps)

    % preallocate: n x 2
    stepSwingStance = zeros(size(legSteps.stepInds,1),2);

    % loop through all steps
    for i = 1:size(legSteps.stepInds, 1)
        % compare absolute translational speed of each half step
        % assign swing/stance correspondingly: faster is swing (-1)
        if (legSteps.stepAbsTransSpd(i,1) >= ...
                legSteps.stepAbsTransSpd(i,2))
            stepSwingStance(i,1) = -1;
            stepSwingStance(i,2) = 1;
        else
            stepSwingStance(i,1) = 1;
            stepSwingStance(i,2) = -1;
        end

    end


    % add swing/stance call to legSteps struct
    legSteps.stepSwingStance = stepSwingStance; 
end