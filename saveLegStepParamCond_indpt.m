% saveLegStepParamCond_indpt.m
%
% Function that saves legStep parameters for all steps that meet the user 
%  specified conditions on bodytraj
% Treats all steps as independent, as long as bodytraj condition(s) are
%  met, unlike saveLegStepParamCond_bouts.m
% Steps are included if at least for half their duration, the bodytraj 
%  condition(s) are met. Step is both swing and stance
% User selects one or more pData files through GUI
% 
% INPUTS:
%   cond - struct of conditions, if multiple conditions, treats it as AND
%       whichParam - cell array (even if 1 element) on which bodytraj field
%           to condition on, one for each condition
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       selLegSteps - all steps, all pData files, that meet conditions
%       selStanceParams - all steps, all pData files, that meet condition,
%           stance only
%       selSwingParams - all steps, all pData files, that meet condition,
%           swing only
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to indices in legSteps) that
%               belong to each pData file, as start and end indices
%       legIDs - legSteps.legIDs
%       cond - same as INPUT
%
% CREATED: 5/8/23 - HHY
%
% UPDATED:
%   5/9/23 - HHY
%
function saveLegStepParamCond_indpt(cond, pDataPath, saveFilePath, ...
    saveFileName)

    % names of all step parameters to save
    stepParamNames = {'stepWhichLeg', 'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY', 'stepAbsDisp', 'stepAbsTransSpd', 'stepT', ...
        'stepBodyFwd', 'stepBodyLat', 'stepBodyYaw'};

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

    % preallocate 
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];

    for i = 1:length(stepParamNames)
        selLegSteps.(stepParamNames{i}) = [];
        selStanceParams.(stepParamNames{i}) = [];
        selSwingParams.(stepParamNames{i}) = [];
    end

    rmvInd = []; % indices of pData files to remove


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
            continue;
        end

        % load legTrack and moveNotMove from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'bodytraj', ...
            'legSteps', 'stanceStepParams', 'swingStepParams');

        % convert legTrack.refPts to legIDs
        legIDs.ind = legTrack.refPts.legInd;
        legIDs.name = legTrack.refPts.legNames;

        % preallocate logical for valid indices of bodytraj
        btValInd = true(size(bodytraj.tZeroed));

        % loop through all conditions
        for j = 1:length(cond.whichParam)
            % check if conditioning on speed, if yes, get absolute value of
            %  corresponding velocity param
            if (contains(cond.whichParam{j}, 'Spd'))
                condWhichParam = replace(cond.whichParam{j}, 'Spd', 'Vel');
                thisCondVar = abs(bodytraj.(condWhichParam));
            % otherwise, just copy over paaram    
            else
                thisCondVar = bodytraj.(cond.whichParam{j});
            end

            % bodytraj indices where this condition is true
            thisValInd = eval(['thisCondVar' cond.cond{j}]);

            % combine with all other conditions
            btValInd = btValInd & thisValInd;
        end

        % filter for moving bouts only
        btValInd(moveNotMove.notMoveInd) = false;

        % preallocate - keep track of indices of steps to include
        inclStepInd = [];

        % loop through all steps of legSteps, check if they're to be
        %  included
        for j = 1:length(legSteps.stepInds)
            % get number of time points for this step, use as approximation
            %  of duration
            numInd = legSteps.stepInds(j,3) - legSteps.stepInds(j,1) + 1;

            numInclInd = sum(btValInd(...
                legSteps.stepInds(j,1):legSteps.stepInds(j,3)));

            % save this step index if more than half of the time during
            %  cond
            if (numInclInd >= numInd/2)
                inclStepInd = [inclStepInd; j];
            end
        end

        % if there are valid steps, save info for this trial
        if (~isempty(inclStepInd))

            % index into selLegSteps, selStanceParams, selSwingParams for start
            %  of steps corresponding to this trial
            trialStartInd = length(selLegSteps.stepWhichLeg) + 1;
    
            % save legStep parameters for this trial
            for j = 1:length(stepParamNames)
                thisLegStep = legSteps.(stepParamNames{j});
                thisLegStepVals = thisLegStep(inclStepInd,:);
                selLegSteps.(stepParamNames{j}) = ...
                    [selLegSteps.(stepParamNames{j}); thisLegStepVals];
    
                thisStanceStep = stanceStepParams.(stepParamNames{j});
                thisStanceVals = thisStanceStep(inclStepInd,:);
                selStanceParams.(stepParamNames{j}) = ...
                    [selStanceParams.(stepParamNames{j}); thisStanceVals];
    
                thisSwingStep = swingStepParams.(stepParamNames{j});
                thisSwingVals = thisSwingStep(inclStepInd,:);
                selSwingParams.(stepParamNames{j}) = ...
                    [selSwingParams.(stepParamNames{j}); thisSwingVals];
            end
    
            % end index for steps for this trial
            trialEndInd = length(selLegSteps.stepWhichLeg);
    
            thisTrialInds = [trialStartInd trialEndInd];
    
            pDataFiles.inds = [pDataFiles.inds; thisTrialInds];
        % otherwise, remove this pData file from list included    
        else
            rmvInd = [rmvInd; i];
        end
    end

    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    % copy over legIDs
    legIDs = legSteps.legIDs;

    % save output file
    fullSavePath = [saveFilePath filesep saveFileName '.mat'];

    save(fullSavePath, 'selLegSteps', 'selStanceParams',...
        'selSwingParams', 'pDataFiles', 'cond', 'legIDs', '-v7.3');
end