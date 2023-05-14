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
%       bouts - struct array, one for each bout 
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
%   5/9/23 - HHY
%
function saveLegStepParamCond_bouts(cond, maxNumSteps, legXYParams, ...
    pDataPath, saveFilePath, saveFileName)

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

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % check if this pData file has legSteps, bodytraj, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'bodytraj')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'bodytraj', ...
            'legSteps', 'stanceStepParams', 'swingStepParams');

        % convert legTrack.refPts to legIDs
        legIDs.ind = legTrack.refPts.legInd;
        legIDs.name = legTrack.refPts.legNames;

        % get yaw velocity peaks for right turns
        [rightPeakInd, rightStartInd, rightEndInd] = findCondYawVelPeaks(...
            bodytraj, cond, moveNotMove, true);
        % get yaw velocity peaks for left turns
        [leftPeakInd, leftStartInd, leftEndInd] = findCondYawVelPeaks(...
            bodytraj, cond, moveNotMove, false);

        
        % get indices for steps aligned to bouts
        % right turns
        [rightStepInd, rightPkSwingStance] = findBoutStepInds(...
            legSteps, rightPeakInd, rightStartInd, rightEndInd, ...
            maxNumSteps, legIDs);
        % left turns
        [leftStepInd, leftPkSwingStance] = findBoutStepInds(...
            legSteps, leftPeakInd, leftStartInd, leftEndInd, ...
            maxNumSteps, legIDs);

        
        % get legStep parameters from aligned step indices

        % get aligned leg X and Y positions from velocity peaks



    end
end