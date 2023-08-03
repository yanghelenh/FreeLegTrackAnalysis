% saveLegStepParamCond_bouts.m
%
% Function that saves legStep parameters for turning bouts, defined by yaw
%  velocity peaks. User can specify conditions in bodytraj that yaw
%  velocity peaks must meet.
% Unlike, saveLegStepParamCond_indpt, extracts bouts that meet bodytraj
%  conditions and accounts for timing of step relative to yaw velocity peak
% Also, saves leg X and Y positions for bouts
% Extracts both right and left turns (flips left turns over midline). Any
%  conditions on yaw or side velocity are considered relative to right
%  turns, and computed for both.
% User selects one or more pData files through GUI
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
%   maxNumSteps - number of steps to each side of peak to consider as part
%       of bout (max bout length is this x2 + 1)
%   legXYParams - struct of parameters for aligning leg X and Y positions
%       maxDuration - time in seconds to consider on each side 
%       interpFrameRate - frame rate to interpolate to
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       selLegSteps - struct of aligned step parameters, where each one is 
%           numSteps x numLegs x 2 x numBouts matrix 
%       selStanceParams - struct of aligned step parameters during stance
%           only, where each one is numSteps x numLegs x numBouts matrix
%       selSwingParams - same as selStanceParams, but for swing
%       pkSwingStance - numLegs x numBouts matrix for whether peak is
%           during swing or stance for each leg and bout
%       stanceParamMeans - struct of means of step parameters during stance
%           only, where each one is numSteps x numLegs matrix 
%       stanceParamStd - as stanceParamMeans, but for standard deviation
%       stanceParamSEM - as stanceParamMeans, but for SEM
%       stanceParamN - as stanceParamMeans, but number of bouts that 
%           contributed to each time point
%       swingParamMeans - as stanceParamMeans, but for swing
%       swingParamStd - as stanceParamStd, but for swing
%       swingParamSEM - as stanceParamSEM, but for swing
%       swingParamN - as stanceParamN, but for swing
%       allLegX - matrix of numTimePts x numLegs x numBouts for aligned 
%           legX positions
%       allLegY - matrix of numTimePts x numLegs x numBouts for aligned
%           legY positions 
%       meanLegX - mean leg X positions, across bouts; numTimePts x 
%           numLegs matrix
%       stdLegX - as meanLegX, but standard deviation
%       SEMLegX - as meanLegX, but SEM
%       meanLegY - mean leg Y positions, across bouts; numTimePts x 
%           numLegs matrix
%       stdLegY - as meanLegY, but standard deviation
%       SEMLegY - as meanLegY, but SEM
%       legT - time points for legX and legY matrices, 0 at yaw vel peak
%       numBouts - total number of bouts (is n for legX and legY, and max n
%           for stepParam (if no NaNs)
%       boutPeakVel - struct for velocity values at each bout peak
%           yaw - vector of length numBouts for peak yaw velocity
%           fwd - vector of length numBouts for peak forward velocity
%           lat - vector of length numBouts for peak lateral velocity
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to bout indices) that
%               belong to each pData file, as cell array of vectors
%       cond - same as INPUT
%       maxNumSteps - same as INPUT
%       legXYParams - same as INPUT
%
% CREATED: 5/8/23 - HHY
%
% UPDATED:
%   5/14/23 - HHY
%   6/6/23 - HHY - update to calculate mean and std appropriately for
%       stepDirections, which is a circular variable. Also, add outlier
%       removal before calculating mean and std
%   6/9/23 - HHY - update to add condition for min and max duration of
%       turning bout (separate from max number of steps)
%   6/14/23 - HHY - fix bug where for stepDirections, time points with no
%       data points returned 0 for the mean and non NaN
%   7/21/23 - HHY - add boutPeakVel to output
%
function saveLegStepParamCond_bouts(cond, maxNumSteps, legXYParams, ...
    pDataPath, saveFilePath, saveFileName)

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY', 'stepAbsDisp', 'stepAbsTransSpd', ...
        'stepBodyFwd', 'stepBodyLat', 'stepBodyYaw'};
    % all the step parameters where values need to be * -1 for left turns
    flipStepParams = {'stepVelY', 'stepAEPY', 'stepPEPY', 'stepBodyLat',...
        'stepBodyYaw', 'stepDirections'};
    % all step parameters that are circular variables - need to use
    %  circular stats - 6/6/23 - HHY
    circStepParams = {'stepDirections'};

    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files', pDataPath, 'MultiSelect', 'on');
    
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
        selLegSteps.(stepParamNames{i}) = [];
        selStanceParams.(stepParamNames{i}) = [];
        selSwingParams.(stepParamNames{i}) = [];
    end

    boutPeakVel.yaw = [];
    boutPeakVel.fwd = [];
    boutPeakVel.lat = [];

    pkSwingStance = [];

    allLegX = [];
    allLegY = [];

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

        % check if this pData file has legSteps, bodytraj, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'bodytraj')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            rmvInd = [rmvInd; i];
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'bodytraj', ...
            'legSteps', 'stanceStepParams', 'swingStepParams');

        % convert legTrack.refPts to legIDs
        legIDs.ind = legTrack.refPts.legInd;
        legIDs.name = legTrack.refPts.legNames;

        % get matching b/w corresponding left and right legs
        rightLegInd = find(contains(legIDs.name, 'R'));
        leftLegInd = find(contains(legIDs.name, 'L'));
        matchedLegInd = zeros(length(rightLegInd),2);

        for j = 1:length(rightLegInd)
            thisLegNum = legIDs.name{rightLegInd(j)}(end);
            thisLeftInd = find(contains(legIDs.name(leftLegInd),thisLegNum));
            matchedLegInd(j,1) = rightLegInd(j);
            matchedLegInd(j,2) = leftLegInd(thisLeftInd);
        end

        % get yaw velocity peaks for right turns
        [rightPeakInd, rightStartInd, rightEndInd] = findCondYawVelPeaks(...
            bodytraj, cond, moveNotMove, true);
        % get yaw velocity peaks for left turns
        [leftPeakInd, leftStartInd, leftEndInd] = findCondYawVelPeaks(...
            bodytraj, cond, moveNotMove, false);


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

        
        % get indices for steps aligned to bouts
        % right turns
        [rightStepInd, rightPkSwingStance] = findBoutStepInds(...
            legSteps, rightPeakInd, rightStartInd, rightEndInd, ...
            maxNumSteps, legIDs);
        % left turns
        [leftStepInd, leftPkSwingStance] = findBoutStepInds(...
            legSteps, leftPeakInd, leftStartInd, leftEndInd, ...
            maxNumSteps, legIDs);

        % peak swing/stance call, merge across legs
        thisPkSwingStance = cat(2, rightPkSwingStance, leftPkSwingStance);
        % add to output matrix
        pkSwingStance = cat(2, pkSwingStance, thisPkSwingStance);



        % get legStep parameters from aligned step indices

        % preallocate
        oneParamValsRight = nan(size(rightStepInd,1), ...
            size(rightStepInd,2), 2, size(rightStepInd,3));
        oneParamStValsRight = nan(size(rightStepInd));
        oneParamSwValsRight = nan(size(rightStepInd));

        oneParamValsLeft = nan(size(leftStepInd,1), ...
            size(leftStepInd,2), 2, size(leftStepInd,3));
        oneParamStValsLeft = nan(size(leftStepInd));
        oneParamSwValsLeft = nan(size(leftStepInd));

        for j = 1:length(stepParamNames)
            thisParam = legSteps.(stepParamNames{j});
            thisStParam = stanceStepParams.(stepParamNames{j});
            thisSwParam = swingStepParams.(stepParamNames{j});
            
            % for right turns
            % loop over all steps of bout
            for k = 1:size(rightStepInd, 1)
                % loop over all legs
                for l = 1:size(rightStepInd, 2)
                    % loop over all bouts
                    for m = 1:size(rightStepInd, 3)
                        thisInd = rightStepInd(k,l,m);
                        % if index is NaN, don't do anything (will keep
                        %  NaN), otherwise, put in appropriate value
                        if ~isnan(thisInd)
                            % legSteps
                            thisVals = thisParam(thisInd,:);
                            oneParamValsRight(k,l,1,m) = thisVals(1);
                            oneParamValsRight(k,l,2,m) = thisVals(2);

                            % stance only
                            oneParamStValsRight(k,l,m) = ...
                                thisStParam(thisInd);
                            % swing only
                            oneParamSwValsRight(k,l,m) = ...
                                thisSwParam(thisInd);
                        end
                    end
                end
            end

            % for left turns
            % loop over all steps of bout
            for k = 1:size(leftStepInd, 1)
                % loop over all legs
                for l = 1:size(leftStepInd, 2)
                    % flip left and right legs
                    % r for row index, c for column index
                    [r, c] = ind2sub(size(matchedLegInd),...
                        find(matchedLegInd == l));
                    % matched ind are across columns, flip columns
                    if (c==1)
                        c = 2;
                    else
                        c = 1;
                    end
                    % new leg index, swapping left and right
                    thisLegInd = matchedLegInd(r,c);

                    % loop over all bouts
                    for m = 1:size(leftStepInd, 3)
                        thisInd = leftStepInd(k,thisLegInd,m);
                        % if index is NaN, don't do anything (will keep
                        %  NaN), otherwise, put in appropriate value
                        if ~isnan(thisInd)
                            % invert values for left turns, for params 
                            %  where it's necessary
                            if any(strcmpi(stepParamNames{j}, flipStepParams))
                                % legSteps
                                thisVals = thisParam(thisInd,:);
                                oneParamValsLeft(k,l,1,m) = ...
                                    thisVals(1) * -1;
                                oneParamValsLeft(k,l,2,m) = ...
                                    thisVals(2) * -1;
    
                                % stance only
                                oneParamStValsLeft(k,l,m) = ...
                                    thisStParam(thisInd) * -1;
                                % swing only
                                oneParamSwValsLeft(k,l,m) = ...
                                    thisSwParam(thisInd) * -1;
                            else
                                % legSteps
                                thisVals = thisParam(thisInd,:);
                                oneParamValsLeft(k,l,1,m) = ...
                                    thisVals(1);
                                oneParamValsLeft(k,l,2,m) = ...
                                    thisVals(2);
    
                                % stance only
                                oneParamStValsLeft(k,l,m) = ...
                                    thisStParam(thisInd);
                                % swing only
                                oneParamSwValsLeft(k,l,m) = ...
                                    thisSwParam(thisInd);
                            end
                        end
                    end
                end
            end

            % concatenate right and left
            oneParamVals = cat(4, oneParamValsRight, oneParamValsLeft);
            oneParamStVals = cat(3, oneParamStValsRight, oneParamStValsLeft);
            oneParamSwVals = cat(3, oneParamSwValsRight, oneParamSwValsLeft);

            % add to output matrices
            selLegSteps.(stepParamNames{j}) = cat(4, ...
                selLegSteps.(stepParamNames{j}), oneParamVals);

            selStanceParams.(stepParamNames{j}) = cat(3, ...
                selStanceParams.(stepParamNames{j}), oneParamStVals);
            selSwingParams.(stepParamNames{j}) = cat(3, ...
                selSwingParams.(stepParamNames{j}), oneParamSwVals);
        end



        % get aligned leg X and Y positions from velocity peaks
        % right turns
        [rightLegX, rightLegY, legT] = getAlignedLegXYfromBouts(legTrack, ...
            legXYParams, rightPeakInd, rightStartInd, rightEndInd);
        % left turns
        [leftLegX, leftLegY, legT] = getAlignedLegXYfromBouts(legTrack, ...
            legXYParams, leftPeakInd, leftStartInd, leftEndInd);

        % flip and invert left turns, but only if there are left turns
        leftLegXFlip = zeros(size(leftLegX));
        leftLegYFlip = zeros(size(leftLegY));
        
        if ~isempty(leftLegX)
            % for left turns, need to swap left and right legs and invert leg Y
            leftLegXFlip(:,matchedLegInd(:,1)) = leftLegX(:,matchedLegInd(:,2));
            leftLegXFlip(:,matchedLegInd(:,2)) = leftLegX(:,matchedLegInd(:,1));
    
            leftLegYFlip(:,matchedLegInd(:,1)) = leftLegY(:,matchedLegInd(:,2));
            leftLegYFlip(:,matchedLegInd(:,2)) = leftLegY(:,matchedLegInd(:,1));
    
            % invert (* -1) leg Y
            leftLegYFlip = leftLegYFlip * -1;
        end

        % concatenate right and flipped/inverted left
        thisLegX = cat(3, rightLegX, leftLegXFlip);
        thisLegY = cat(3, rightLegY, leftLegYFlip);

        % add to total output matrix
        % if it's the first trial to add data
        if isempty(allLegX)
            allLegX = thisLegX;
            allLegY = thisLegY;
        % all other trials are appended to end    
        else
            allLegX = cat(3, allLegX, thisLegX);
            allLegY = cat(3, allLegY, thisLegY);
        end

    end

    % compute means, std dev, SEM, n on selStanceParams, selSwingParams
    % preallocate
    for i = 1:length(stepParamNames)
        stanceParamMeans.(stepParamNames{i}) = nan(...
            size(selStanceParams.(stepParamNames{i}),1),...
            size(selStanceParams.(stepParamNames{i}),2));
        stanceParamStd.(stepParamNames{i}) = nan(...
            size(selStanceParams.(stepParamNames{i}),1),...
            size(selStanceParams.(stepParamNames{i}),2));
        stanceParamSEM.(stepParamNames{i}) = nan(...
            size(selStanceParams.(stepParamNames{i}),1),...
            size(selStanceParams.(stepParamNames{i}),2));
        stanceParamN.(stepParamNames{i}) = nan(...
            size(selStanceParams.(stepParamNames{i}),1),...
            size(selStanceParams.(stepParamNames{i}),2));

        swingParamMeans.(stepParamNames{i}) = nan(...
            size(selSwingParams.(stepParamNames{i}),1),...
            size(selSwingParams.(stepParamNames{i}),2));
        swingParamStd.(stepParamNames{i}) = nan(...
            size(selSwingParams.(stepParamNames{i}),1),...
            size(selSwingParams.(stepParamNames{i}),2));
        swingParamSEM.(stepParamNames{i}) = nan(...
            size(selSwingParams.(stepParamNames{i}),1),...
            size(selSwingParams.(stepParamNames{i}),2));
        swingParamN.(stepParamNames{i}) = nan(...
            size(selSwingParams.(stepParamNames{i}),1),...
            size(selSwingParams.(stepParamNames{i}),2));
    end
    
    % loop through all params
    for i = 1:length(stepParamNames)
        % for stance
        thisStanceParam = selStanceParams.(stepParamNames{i});
        % loop through all time points
        for j = 1:size(thisStanceParam,1)
            % loop through all legs
            for k = 1:size(thisStanceParam,2)
                thisTandLegs = thisStanceParam(j,k,:);

                % remove NaNs
                thisTandLegs(isnan(thisTandLegs)) = [];

                % 6/6/23 - add outlier removal before calculating mean and
                %  std
                % get mean and std dev, account for if circular 
                if(any(strcmpi(stepParamNames{i}, circStepParams)))
                    % convert this parameter to radians
                    thisTandLegs = deg2rad(thisTandLegs);
                    % get circular mean
                    thisMean = rad2deg(circ_mean(rmoutliers(squeeze(thisTandLegs))));
                    % get circular std
                    thisStd = rad2deg(circ_std(rmoutliers(squeeze(thisTandLegs))));
                % if not circular, compute regular mean and std
                else
                    thisMean = mean(rmoutliers(squeeze(thisTandLegs)));
                    thisStd = std(rmoutliers(squeeze(thisTandLegs)));
                end

                % get SEM, n for this time point and legs
                thisN = length(thisTandLegs);
                thisSEM = thisStd / sqrt(thisN);

                % add to output matrices if there are any data points
                if (thisN > 0)
                    stanceParamMeans.(stepParamNames{i})(j,k) = thisMean;
                    stanceParamStd.(stepParamNames{i})(j,k) = thisStd;
                    stanceParamSEM.(stepParamNames{i})(j,k) = thisSEM;
                    stanceParamN.(stepParamNames{i})(j,k) = thisN;
                else
                    stanceParamMeans.(stepParamNames{i})(j,k) = nan;
                    stanceParamStd.(stepParamNames{i})(j,k) = nan;
                    stanceParamSEM.(stepParamNames{i})(j,k) = nan;
                    stanceParamN.(stepParamNames{i})(j,k) = thisN;
                end
            end
        end

        % for swing
        thisSwingParam = selSwingParams.(stepParamNames{i});
        % loop through all time points
        for j = 1:size(thisSwingParam,1)
            % loop through all legs
            for k = 1:size(thisSwingParam,2)
                thisTandLegs = thisSwingParam(j,k,:);

                % remove NaNs
                thisTandLegs(isnan(thisTandLegs)) = [];

                % get mean and std dev, account for if circular 
                if(any(strcmpi(stepParamNames{i}, circStepParams)))
                    % convert this parameter to radians
                    thisTandLegs = deg2rad(thisTandLegs);
                    % get circular mean
                    thisMean = rad2deg(circ_mean(rmoutliers(squeeze(thisTandLegs))));
                    % get circular std
                    thisStd = rad2deg(circ_std(rmoutliers(squeeze(thisTandLegs))));
                % if not circular, compute regular mean and std
                else
                    thisMean = mean(rmoutliers(squeeze(thisTandLegs)));
                    thisStd = std(rmoutliers(squeeze(thisTandLegs)));
                end

                % get SEM, n for this time point and legs
                thisN = length(thisTandLegs);
                thisSEM = thisStd / sqrt(thisN);

                % add to output matrices
                swingParamMeans.(stepParamNames{i})(j,k) = thisMean;
                swingParamStd.(stepParamNames{i})(j,k) = thisStd;
                swingParamSEM.(stepParamNames{i})(j,k) = thisSEM;
                swingParamN.(stepParamNames{i})(j,k) = thisN;
            end
        end
    end

    % compute means, std dev, SEM, n on allLegX and allLegY
    meanLegX = mean(allLegX, 3);
    stdLegX = std(allLegX, [], 3);
    SEMLegX = stdLegX / sqrt(countNumBouts);

    meanLegY = mean(allLegY, 3);
    stdLegY = std(allLegY, [], 3);
    SEMLegY = stdLegY / sqrt(countNumBouts); 


    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    numBouts = countNumBouts;

    % save output file
    fullSavePath = [saveFilePath filesep saveFileName '.mat'];

    save(fullSavePath, 'selLegSteps', 'selStanceParams', ...
        'selSwingParams', 'pkSwingStance', 'stanceParamMeans', ...
        'stanceParamStd', 'stanceParamSEM', 'stanceParamN', ...
        'swingParamMeans', 'swingParamStd', 'swingParamSEM', ...
        'swingParamN', 'allLegX', 'allLegY', 'meanLegX', 'stdLegX', ...
        'SEMLegX', 'meanLegY', 'stdLegY', 'SEMLegY', 'legT', 'numBouts', ...
        'boutPeakVel', 'pDataFiles', 'cond', 'maxNumSteps', ...
        'legXYParams', '-v7.3');
end