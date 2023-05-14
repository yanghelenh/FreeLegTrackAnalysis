% extractStepsFromPData.m
%
% Function that computes leg step parameters from raw leg tracking data 
%  (output of loadFreeWalkTrk2PData())
% User selects one or more pData files through GUI
% The pData file must contain the moveNotMove struct (output of
%  extractMoveNotMoveFromPData() and the legTrack struct 
%  (output of loadFreeWalkTrk2PData())
% Saves back into same pData file new structs legSteps, stanceStepParams,
%  and swingStepParams
%
% INPUTS:
%   pDataPath - full path to pData folder
%   legRevParams - struct of parameters to use for getting leg position reversals
%       if empty, use default parameters
%
% OUTPUTS: none, but save legSteps, stanceStepParams, and swingStepParams
%   back into pData file
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%
function extractStepsFromPData(pDataPath, legRevParams)

    % default parameters for extracting peaks
    if isempty(legRevParams)
        legRevParams.minProm = 0.05; % MinPeakProminence of findpeaks
        legRevParams.minDist = 6; % MinPeakDistance of findpeaks
    end


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
    
        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legTrack struct, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legTrack')))
            fprintf('%s does not contain the legTrack struct. Ending \n', ...
                pDataName);
            continue;
        elseif (~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            fprintf('%s does not contain the moveNotMove struct. Ending \n', ...
                pDataName);
            continue;
        end

        % load legTrack and moveNotMove from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'bodytraj');

        % convert legTrack.refPts to legIDs
        legIDs.ind = legTrack.refPts.legInd;
        legIDs.name = legTrack.refPts.legNames;


        % get indices of max and min for 6 legs
        [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
            getLegReversals(legTrack, moveNotMove, legRevParams, legIDs);

        % add these parameters to legSteps struct
        legSteps.maxIndsAll = maxIndsAll;
        legSteps.minIndsAll = minIndsAll;
        legSteps.maxWhichLeg = maxWhichLeg;
        legSteps.minWhichLeg = minWhichLeg;
        legSteps.legIDs = legIDs;

        % compute steps
        [stepInds, stepWhichLeg] = minMaxPosToSteps(legSteps.maxIndsAll, ...
            legSteps.minIndsAll, legSteps.maxWhichLeg, legSteps.minWhichLeg,...
            moveNotMove.notMoveBout, moveNotMove.moveBout);
        
        % save into legSteps struct
        legSteps.stepInds = stepInds;
        legSteps.stepWhichLeg = stepWhichLeg;
        
        
        % compute step parameters
        legSteps = computeFreeWalkStepParameters(legSteps, legTrack, ...
            bodytraj);
        
        % get swing/stance - based on movement in arena reference frame
        legSteps = callFreeWalkSwingStanceSteps(legSteps);
        
        % step parameters during swing/stance
        [stanceStepParams, swingStepParams] = getStepParamsSwingStance(...
            legSteps);
    
        % update pData file
        save(pDataFullPath, 'legSteps', 'stanceStepParams', ...
            'swingStepParams', '-append');

    end
end