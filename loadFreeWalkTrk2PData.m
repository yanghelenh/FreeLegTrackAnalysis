% loadFreeWalkTrk2PData.m
%
% Function that finds the .trk file corresponding to the selected pData
%  file(s), preprocesses the .trk data, and saves it into the pData file
% Assumes that .trk file name matches pData name in trial-#-date-time
% .trk file as standard output from APT
%
% Modified version of loadTrk2PData (from tethered fly tracking)
%
% NOTE: 
%
% INPUTS:
%   pDataPath - full path to folder containing pData files
%   trkFilePath - full path to all .trk files
%   trkSuffix - "v#_mdn_date.trk" - version number and model type and date
%       of last training .trk suffix for .trk files; specify in case there 
%       are multiple .trk files for a given trial
%   refPts - struct of parameters defining reference tracked points
%       (optional, pass [] to use defaults)
%   smoParams - struct of parameters for smoothing leg postion and 
%       velocity (optional, pass [] to use defaults)
%   
% OUTPUTS:
%   none, but saves processed .trk data into pData file
%
% CREATED: 4/8/23 - HHY
%
% UPDATED:
%   4/9/23 - HHY
%
function loadFreeWalkTrk2PData(pDataPath, trkFilePath, trkSuffix, ...
    refPts, smoParams)

    % SOME DEFAULT CONSTANTS %
    % if refPts not specified, default params
    if (isempty(refPts))
        % indicies for specific body parts
        refPts.legInd = 1:6;
        refPts.legNames = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'};
        refPts.headPtInd = 7;
        refPts.thxRPtInd = 8;
        refPts.thxLPtInd = 9;
        refPts.thxMidInd = 10;
        refPts.abdPtInd = 11;
    end

    % constants
    % conversion between camera pixels and meters
    imgParams.PX_PER_METER = 37023.1016957; 
    % number of pixels in x dimension of video frame
    imgParams.xPx = 658;
    % number of pixels in y dimension of video frame
    imgParams.yPx = 494;

    % if smoParams not specified, default params
    if (isempty(smoParams))
        % parameters for smoothing for determining leg velocities
        smoParams.padLen = 50; % pad length in samples
        smoParams.sigma = 5; % in samples
        % parameters for removing outliers
        smoParams.medFiltDeg = 5; % degree of median filter
        % percent of (max-min) is acceptable as deviation
        smoParams.percntDev = 50; 
        smoParams.maxMinWin = 500; % size of window to compute max and min
    end


    % check that trkFilePath does not end with filesep
    if (strcmp(trkFilePath(end),filesep))
        % remove if last character is filesep
        trkFilePath = trkFilePath(1:(end-1));
    end

    % get names of all .trk files in folder specified by trkFilePath with
    %  the trkSuffix
    trkFilesSearch = [trkFilePath filesep '*' trkSuffix];
    trkFilesStrct = dir(trkFilesSearch);

    % convert files struct into cell array of names
    trkFilesNames = cell(size(trkFilesStrct));

    for i = 1:length(trkFilesStrct)
        trkFilesNames{i} = trkFilesStrct(i).name;
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

        % extract file name without .mat suffix
        % assumes format of trial-#-date-time.mat
        trialFullName = pDataName(1:(end-4));

        % full path for this pData file
        pDataFullPath = [pDataDirPath pDataName];
        
        % load bodytraj struct
        load(pDataFullPath, 'bodytraj');
        
        % logical of whether trial full name matches names in trk files
        %  list
        trkLog = contains(trkFilesNames, trialFullName);

        % only proceed if there is a matching .trk file
        if (any(trkLog))
            % get .trk file name
            % assumes only 1 matching .trk file (which should be the
            % case)
            matchedTrkFileName = trkFilesNames{trkLog};

            % .trk full path
            trkFileFullPath = [trkFilePath filesep matchedTrkFileName];
            

            % preprocess .trk data
            legTrack = preprocessFreeWalkLegTrack(trkFileFullPath, ...
                bodytraj.tZeroed, bodytraj.cam.startFr, ...
                bodytraj.cam.endFr, bodytraj.cnc.x, bodytraj.cnc.y, ...
                refPts, smoParams, imgParams);

            % save legTrack struct into pData
            save(pDataFullPath, 'legTrack', '-append'); 

            % display
            fprintf('Saved legTrack for %s!\n', pDataName);
        else
            fprintf('%s does not have a matching .trk file\n', ...
                pDataName);
        end
    end
end