% findBoutStepInds.m
%
% Helper function for saveLegStepParamCond_bouts() that takes in bout peak,
%  start, and end indices (output of findCondYawVelPeaks() helper function)
%  as well as step index info, and returns aligned steps, as indices into
%  legStep matrices.
% Also, returns whether bout peak is during stance or swing for each leg
%
% INPUTS:
%   legSteps - struct of legStep info
%   peakInd - indices of frames of bout peaks
%   boutStartInd - indices of bout starts
%   boutEndInd - indices of bout ends
%   maxNumSteps - number of steps to each side of peak to consider as part
%       of bout (max bout length is this x2 + 1)
%   legIDs - struct defining indices corresponding to each leg
%
% OUTPUTS:
%   boutStepInd - maxNumSteps*2+1 x numLegs x numBouts size matrix defining
%       aligned bouts as indices into legSteps matrices
%   boutPkSwingStance - numLegs x numBouts size matrix indicating whether
%       the peak of each bout is during swing or stance for each leg (-1
%       for swing, +1 for stance)
%
% CREATED: 5/11/23 - HHY
%
% UPDATED:
%   5/11/23 - HHY
%
function [boutStepInd, boutPkSwingStance] = findBoutStepInds(legSteps, ...
    peakInd, boutStartInd, boutEndInd, maxNumSteps, legIDs)

    % some params
    numLegs = length(legIDs.ind);
    numBouts = length(peakInd);

    % preallocate
    boutStepInd = nan(maxNumSteps * 2 + 1, numLegs, numBouts);
    boutPkSwingStance = nan(numLegs, numBouts);

    % loop through all peaks
    for i = 1:length(peakInd)
        thisPeakInd = peakInd(i); 

        % loop through all legs
        for j = 1:numLegs
            % logical for steps for this leg only
            thisLegLog = legSteps.stepWhichLeg == legIDs.ind(j);

            % get index of step during which peak occurs
            % if peak is on step boundary, takes later step (based on which
            %  one is equal)
            % also ensure it's for this leg only
            stepAtPeak = find((legSteps.stepInds(:,1)<=thisPeakInd) & ...
                (legSteps.stepInds(:,3)>thisPeakInd) & thisLegLog);

            % check that a step is present
            % if a step is not present, skip to next leg
            if (isempty(stepAtPeak))
                continue;
            end

            % determine whether peak falls during step's swing or stance
            % if the peak is during the first half
            if (thisPeakInd < legSteps.stepInds(stepAtPeak, 2))
                thisSwSt = legSteps.stepSwingStance(stepAtPeak, 1);
            % otherwise, in second half    
            else
                thisSwSt = legSteps.stepSwingStance(stepAtPeak, 2);
            end
            boutPkSwingStance(j,i) = thisSwSt;

            % place step at peak into output matrix
            boutStepInd(maxNumSteps + 1, j, i) = stepAtPeak;

            % get steps prior to peak step, place in output matrix
            for k = 1:maxNumSteps
                % adjacent steps should have adjacent indices in stepInds
                %  and end of one step should equal start of another in
                %  frame indices
                adjStepInd = stepAtPeak - k;

                % if the step index is before trial start, end update
                if (adjStepInd < 1)
                    break;
                % check that adjacent step still belongs to same leg    
                elseif (legSteps.stepWhichLeg(adjStepInd)~=legIDs.ind(j))
                    break;
                % check that end frame of this step matches start of next one
                % if not, end update
                elseif (legSteps.stepInds(adjStepInd,3)~=...
                        legSteps.stepInds(adjStepInd+1,1))
                    break;
                % check that the start of this step happens after the start 
                %  the bout
                elseif (legSteps.stepInds(adjStepInd,1) < ...
                        boutStartInd(i))
                    break;
                % if valid, add to matrix    
                else
                    boutStepInd(maxNumSteps + 1 - k,j,i) = adjStepInd;
                end
            end

            % get steps after peak step, place in output matrix
            for k = 1:maxNumSteps
                % adjacent steps should have adjacent indices in stepInds
                %  and end of one step should equal start of another in
                %  frame indices
                adjStepInd = stepAtPeak + k;

                % if the step index is after trial end, end update
                if (adjStepInd > length(legSteps.stepWhichLeg))
                    break;
                % check that adjacent step still belongs to same leg 
                elseif (legSteps.stepWhichLeg(adjStepInd)~=legIDs.ind(j))
                    break;
                % check that end of last frame matches start of this one
                % otherwise, end update
                elseif (legSteps.stepInds(adjStepInd,1) ~= ...
                        legSteps.stepInds(adjStepInd - 1,3))
                    break;
                % check that the end of this step happens before the end of
                %  the bout
                elseif (legSteps.stepInds(adjStepInd,3) > ...
                        boutEndInd(i))
                    break;
                % if valid, add to matrix
                else
                    boutStepInd(maxNumSteps + 1 + k, j, i) = adjStepInd;
                end
            end

        end
    end
end