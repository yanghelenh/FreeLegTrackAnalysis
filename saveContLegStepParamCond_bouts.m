% saveContLegStepParamCond_bouts.m
%
% Function that saves continuous leg step parameters, aligned to yaw
%  velocity peak. User specifies conditions that bodytraj parameters must
%  meet. 
% Adaptation of saveLegStepParamCond_bouts(), except it uses legStepsCont
%  instead of legSteps step parameters. I.e. it uses the continuous
%  estimate of the legStep parameters
%
% INPUTS:
%   cond - struct of conditions for yaw velocity peak, if multiple 
%     conditions, treats it as AND
%       whichParam - cell array (even if 1 element) on which bodytraj field
%           to condition on, one for each condition
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%       turnDur - 2 element vector [minTurnDuration maxTurnDuration] to
%           specify the min and max duration of the turning bout for it to
%           be included
%       minYawThresh - minimum yaw velocity to define start and end of bout
%   fwdVelCond - struct of conditions on forward velocity to apply to
%     turning bout
%       initVel - 2 element vector defining [min max] range of acceptable
%           initial forward velocities (at bout start)
%       change - 2 element vector defining [min max] range of acceptable
%           changes in forward velocity, peak - start 
%   validDur - time in seconds on each side of peak to consider
%   interpFramerate - frame rate to interpolate output continuous leg
%       parameters
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       selLegStepsCont - struct of aligned continuous step parameters,
%           where each one is numTimePts x numLegs x numBouts matrix
%       legStepsContMean - struct of means over all bouts for each step 
%           parameter, as numTimePts x numLegs matrix for each
%       legStepsContStd - struct of std dev over all bouts for each step
%           parameter, as numTimePts x numLegs matrix for each
%       legStepsContSEM - struct of SEMs over all bouts for each step
%           parameter, as numTimePts x numLegs matrix for each
%       legT - time points for matrices, 0 at yaw vel peak
%       numBouts - total number of bouts
%       boutPeakVel - struct for velocity values at each bout peak
%           yaw - vector of length numBouts for peak yaw velocity
%           fwd - vector of length numBouts for peak forward velocity
%           lat - vector of length numBouts for peak lateral velocity
%       boutStartVel - struct for velocity values at each bout start
%           yaw - vector of length numBouts for start yaw velocity
%           fwd - vector of length numBouts for start forward velocity
%           lat - vector of length numBouts for start lateral velocity
%       boutEndVel - struct for velocity values at each bout end
%           yaw - vector of length numBouts for end yaw velocity
%           fwd - vector of length numBouts for end forward velocity
%           lat - vector of length numBouts for end lateral velocity
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to bout indices) that
%               belong to each pData file, as cell array of vectors
%       cond - same as INPUT
%       fwdVelCond - same as INPUT
%
% CREATED: 9/5/23 - HHY
%
% UPDATED:
%   9/5/23 - HHY
%
function saveContLegStepParamCond_bouts(cond, fwdVelCond, validDur, ...
    interpFramerate, pDataFNames, pDataPath, saveFilePath, saveFileName)

    % names of all step parameters to save
    stepParamNames = {'AEPX', 'PEPX', 'AEPY', 'PEPY', 'stepLengthX', ...
        'stepLengthY', 'stepLength', 'stepDirection'};
    % all the step parameters where values need to be * -1 for left turns
    flipStepParams = {'AEPY', 'PEPY', 'stepLengthY', 'stepDirection'};
    % all step parameters that are circular variables - need to use
    %  circular stats
    circStepParams = {'stepDirection'};

    % prompt user to select pData files
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files', pDataPath, 'MultiSelect', 'on');
    else
        pDataDirPath = pDataPath;
    end
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate outputs
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];

    for i = 1:length(stepParamNames)
        selLegStepsCont.(stepParamNames{i}) = [];
    end

    boutPeakVel.yaw = [];
    boutPeakVel.fwd = [];
    boutPeakVel.lat = [];

    boutStartVel.yaw = [];
    boutStartVel.fwd = [];
    boutStartVel.lat = [];

    boutEndVel.yaw = [];
    boutEndVel.fwd = [];
    boutEndVel.lat = [];


    rmvInd = []; % indices of pData files to remove
    countNumBouts = 0; % counter for number of bouts

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legStepsCont, bodytraj, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'bodytraj')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            rmvInd = [rmvInd; i];
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'bodytraj', ...
            'legSteps', 'legStepsCont');

        % convert legTrack.refPts to legIDs
        legIDs = legStepsCont.legIDs;

        % get matching b/w corresponding left and right legs
        rightLegInd = find(contains(legIDs.names, 'R'));
        leftLegInd = find(contains(legIDs.names, 'L'));
        matchedLegInd = zeros(length(rightLegInd),2);

        for j = 1:length(rightLegInd)
            thisLegNum = legIDs.names{rightLegInd(j)}(end);
            thisLeftInd = find(contains(legIDs.names(leftLegInd),thisLegNum));
            matchedLegInd(j,1) = rightLegInd(j);
            matchedLegInd(j,2) = leftLegInd(thisLeftInd);
        end

        % get yaw velocity peaks for right turns
        [rightPeakInd, rightStartInd, rightEndInd] = findCondYawVelPeaks(...
            bodytraj, cond, fwdVelCond, moveNotMove, true);
        % get yaw velocity peaks for left turns
        [leftPeakInd, leftStartInd, leftEndInd] = findCondYawVelPeaks(...
            bodytraj, cond, fwdVelCond, moveNotMove, false);


        % check if this pData file contributes any turns
        % if not, skip and move to next pData file
        if (isempty(rightPeakInd) && isempty(leftPeakInd))
            rmvInd = [rmvInd;i];
            continue;
        end

        % number of bouts for this trial
        thisNumBouts = length(rightPeakInd) + length(leftPeakInd);

        % get indices for bouts for this trial
        thisTrialStartInd = 1 + countNumBouts;
        thisTrialEndInd = thisTrialStartInd + thisNumBouts - 1;

        thisTrialInds = [thisTrialStartInd thisTrialEndInd];
    
        pDataFiles.inds = [pDataFiles.inds; thisTrialInds];

        % update counter
        countNumBouts = countNumBouts + thisNumBouts;

        % get smoothed bodytraj values for peaks
        boutPeakVel.yaw = [boutPeakVel.yaw; ...
            bodytraj.angVelSmoS(rightPeakInd)];
        boutPeakVel.yaw = [boutPeakVel.yaw; ...
            bodytraj.angVelSmoS(leftPeakInd)];
        boutPeakVel.fwd = [boutPeakVel.fwd; ...
            bodytraj.fwdVelSmoS(rightPeakInd)];
        boutPeakVel.fwd = [boutPeakVel.fwd; ...
            bodytraj.fwdVelSmoS(leftPeakInd)];
        boutPeakVel.lat = [boutPeakVel.lat; ...
            bodytraj.latVelSmoS(rightPeakInd)];
        boutPeakVel.lat = [boutPeakVel.lat; ...
            bodytraj.latVelSmoS(leftPeakInd)];

        % get smoothed bodytraj values for bout starts
        boutStartVel.yaw = [boutStartVel.yaw; ...
            bodytraj.angVelSmoS(rightStartInd)];
        boutStartVel.yaw = [boutStartVel.yaw; ...
            bodytraj.angVelSmoS(leftStartInd)];
        boutStartVel.fwd = [boutStartVel.fwd; ...
            bodytraj.fwdVelSmoS(rightStartInd)];
        boutStartVel.fwd = [boutStartVel.fwd; ...
            bodytraj.fwdVelSmoS(leftStartInd)];
        boutStartVel.lat = [boutStartVel.lat; ...
            bodytraj.latVelSmoS(rightStartInd)];
        boutStartVel.lat = [boutStartVel.lat; ...
            bodytraj.latVelSmoS(leftStartInd)];

        % get smoothed bodytraj values for bout ends
        boutEndVel.yaw = [boutEndVel.yaw; ...
            bodytraj.angVelSmoS(rightEndInd)];
        boutEndVel.yaw = [boutEndVel.yaw; ...
            bodytraj.angVelSmoS(leftEndInd)];
        boutEndVel.fwd = [boutEndVel.fwd; ...
            bodytraj.fwdVelSmoS(rightEndInd)];
        boutEndVel.fwd = [boutEndVel.fwd; ...
            bodytraj.fwdVelSmoS(leftEndInd)];
        boutEndVel.lat = [boutEndVel.lat; ...
            bodytraj.latVelSmoS(rightEndInd)];
        boutEndVel.lat = [boutEndVel.lat; ...
            bodytraj.latVelSmoS(leftEndInd)];

        

        % loop through all step parameters, get bout vals
        for j = 1:length(stepParamNames)
            % get values for this step parameter, right turn bouts
            [oneParamValsRight, legT] = getAlignedLegStepsContfromBouts(...
                stepParamNames{j}, legStepsCont, validDur, ...
                interpFramerate, rightPeakInd, false, matchedLegInd, false);

            % check whether this is a parameter where the value needs to be
            %  inverted for left turns
            if any(strcmpi(stepParamNames{j}, flipStepParams))
                thisInv = true;
            else
                thisInv = false;
            end

            % get values for this step parameter, left turn bouts
            [oneParamValsLeft, ~] = getAlignedLegStepsContfromBouts(...
                stepParamNames{j}, legStepsCont, validDur, ...
                interpFramerate, leftPeakInd, true, matchedLegInd, thisInv);

             % concatenate right and left
            oneParamVals = cat(3, oneParamValsRight, oneParamValsLeft);

            % add to output matrix
            selLegStepsCont.(stepParamNames{j}) = cat(3, ...
                selLegStepsCont.(stepParamNames{j}), oneParamVals);
        end
    end

    % compute means, std dev, SEM on selLegStepsCont
    % preallocate
    for i = 1:length(stepParamNames)
        legStepsContMean.(stepParamNames{i}) = nan(...
            size(selLegStepsCont.(stepParamNames{i}),1),...
            size(selLegStepsCont.(stepParamNames{i}),2));
        legStepsContStd.(stepParamNames{i}) = nan(...
            size(selLegStepsCont.(stepParamNames{i}),1),...
            size(selLegStepsCont.(stepParamNames{i}),2));
        legStepsContSEM.(stepParamNames{i}) = nan(...
            size(selLegStepsCont.(stepParamNames{i}),1),...
            size(selLegStepsCont.(stepParamNames{i}),2));
    end
    
    % loop through all params
    for i = 1:length(stepParamNames)
        thisParamVal = selLegStepsCont.(stepParamNames{i});
        % loop through all time points
        for j = 1:size(thisParamVal,1)
            % loop through all legs
            for k = 1:size(thisParamVal,2)
                thisTandLegs = thisParamVal(j,k,:);

                % get mean and std dev, account for if circular 
                if(any(strcmpi(stepParamNames{i}, circStepParams)))
                    % convert this parameter to radians
                    thisTandLegs = deg2rad(thisTandLegs);
                    % get circular mean
                    thisMean = rad2deg(circ_mean(squeeze(thisTandLegs)));
                    % get circular std
                    thisStd = rad2deg(circ_std(squeeze(thisTandLegs)));
                % if not circular, compute regular mean and std
                else
                    thisMean = mean(squeeze(thisTandLegs));
                    thisStd = std(squeeze(thisTandLegs));
                end

                % get SEM, n for this time point and legs
                thisN = length(thisTandLegs);
                thisSEM = thisStd / sqrt(thisN);

                % add to output matrices
                legStepsContMean.(stepParamNames{i})(j,k) = thisMean;
                legStepsContStd.(stepParamNames{i})(j,k) = thisStd;
                legStepsContSEM.(stepParamNames{i})(j,k) = thisSEM;
            end
        end
    end

    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    numBouts = countNumBouts;

    % save output file
    fullSavePath = [saveFilePath filesep saveFileName '.mat'];

    save(fullSavePath, 'selLegStepsCont', 'legStepsContMean', ...
        'legStepsContStd', 'legStepsContSEM', 'legT', 'numBouts', ...
        'boutPeakVel', 'boutStartVel', 'boutEndVel', 'pDataFiles', ...
        'cond', 'fwdVelCond', '-v7.3');
end
