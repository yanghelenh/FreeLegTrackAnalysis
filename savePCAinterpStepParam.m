% savePCAinterpStepParam.m
%
% Function that takes in two sets of step parameters, performs PCA on each,
%  and compares the two against each other.
% Unlike pcaCompStepParams(), this function does not operate on turn bouts.
%  Instead, it works on step parameters interpolated to the specified frame
%  rate.
% Select pData files through GUI
% Saves output file with value of PC1 for each set for all time points.
%  Also saves values of all other step parameters and bodytraj velocity 
%  values at these same time points.
%
% INPUTS:
%   set1 - struct of first set of step params
%     name - string for name of this set
%     params - cell array of step param names, 1 per variable
%     legs - cell array of leg names (R1-3, L1-3), 1 per variable
%     whichPhase - which phase ('swing' or 'stance')
%   set2 - struct of second set of step params
%     name - string for name of this set
%     params - cell array of step param names, 1 per variable
%     legs - cell array of leg names (R1-3, L1-3), 1 per variable
%     whichPhase - which phase ('swing' or 'stance')
%   interpFrameRate - sampling rate to interpolate step parameters to, in
%       Hz
%   pDataPath - full path to pData directory
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   none, but saves output file with following variables
%       set1 - struct of first set of step parameters, from input, with
%           additional fields:
%         coeff - PCA coefficients
%         latent - PCA latent
%         tsquared - PCA tsquared
%         explained - PCA variance explained
%         mu - mean for each dimension
%       set2 - struct of second set of step parameters, from input, with
%           additional fields as set 1
%       set1PCVars - numSampPts x numSet1Vars matrix that PCA is run on
%       set2PCVars - numSampPts x numSet2Vars matrix that PCA is run on
%       set1Score - PC values for set1 
%       set2Score - PC values for set2
%       matchStanceStepParams - step parameter values during stance,
%         at same time points as PC scores
%       matchSwingStepParams - step parameter values during swing, at 
%         same time points as PC scores
%       matchBodytrajParams - body trajectory values at same time points as
%         PC scores
%       legIDs - leg info, copied from legTrack
%         ind - leg indices
%         name - leg names
%       pDataFiles - info on pData files that contribute to the output
%         names - file names
%         inds - start and end indices for points for each file
%       interpFrameRate - copied from INPUTS
%
% CREATED: 7/30/23 - HHY
%
% UPDATED:
%   7/30/23 - HHY
%   8/2/23 - HHY - update comments
%
function savePCAinterpStepParam(set1, set2, interpFrameRate, pDataPath, ...
    saveFileName, saveFileDir)

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY'};

    bodytrajParamNames = {'angVelSmoM', 'fwdVelSmoM', 'latVelSmoM'};

    % number of step parameter variables
    set1NumVars = length(set1.params);
    set2NumVars = length(set2.params);

    % interpolation ifi
    ifi = 1/interpFrameRate;

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
    rmvPDataInd = [];
    currNumPts = 0;

    % matrices on which to run PCA; each column is variable, each row is a
    %  time point; pooled across pData files
    set1PCVars = [];
    set2PCVars = [];

    % for step parameters at matched data points to PC scores
    for i = 1:length(stepParamNames)
        matchStanceStepParams.(stepParamNames{i}) = [];
        matchSwingStepParams.(stepParamNames{i}) = [];
    end

    % for bodytraj parameters at matched data points
    for i = 1:length(bodytrajParamNames)
        matchBodytrajParams.(bodytrajParamNames{i}) = [];
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

        % check if this pData file has legSteps, bodytraj, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'bodytraj')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            rmvPDataInd = [rmvPDataInd; i];
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'bodytraj', ...
            'legSteps', 'stanceStepParams', 'swingStepParams');

        % convert legTrack.refPts to legIDs
        legIDs.ind = legTrack.refPts.legInd;
        legIDs.name = legTrack.refPts.legNames;

        % get sample times to interpolate to
        startTime = legTrack.t(1);
        endTime = legTrack.t(end);

        sampTime = startTime:ifi:endTime;
        sampTime = sampTime';

        % check that there are at least 2 steps for each leg for this pData
        %  file
        tooFewSteps = false;
        for j = 1:length(legIDs.ind)
            if (sum(legSteps.stepWhichLeg == legIDs.ind(j)) < 2)
                tooFewSteps = true;
            end
        end

        % if not, skip to next pData file
        if (tooFewSteps)
            rmvPDataInd = [rmvPDataInd; i];
            continue;
        end

        % preallocate matrix of var vals for PCA
        pcaVarsMat = zeros(length(sampTime), set1NumVars + set2NumVars);
        
        % loop through all set 1 variables
        for j = 1:set1NumVars
            pcaVarsMat(:,j) = stepParam2VectorSpline(legSteps, ...
                legTrack.t, sampTime, ...
                moveNotMove, set1.legs{j}, set1.params{j}, ...
                set1.whichPhase{j});
        end

        % loop through all set 2 variables
        for j = 1:set2NumVars
            pcaVarsMat(:,j+set1NumVars) = stepParam2VectorSpline(...
                legSteps, legTrack.t, sampTime, moveNotMove, ...
                set2.legs{j}, set2.params{j}, set2.whichPhase{j});
        end

        % get interpolated values for all step parameters
        for j = 1:length(stepParamNames)
            % reinialize for each step parameter
            thisStanceStepVal = zeros(length(sampTime), length(legIDs.ind));
            thisSwingStepVal = zeros(length(sampTime), length(legIDs.ind));

            % loop through all legs
            for k = 1:length(legIDs.ind)
                thisStanceStepVal (:,k) = ...
                stepParam2VectorSpline(legSteps, legTrack.t, sampTime, ...
                moveNotMove, legIDs.name{k}, stepParamNames{j}, ...
                'stance');

                thisSwingStepVal (:,k) = ...
                stepParam2VectorSpline(legSteps, legTrack.t, sampTime, ...
                moveNotMove, legIDs.name{k}, stepParamNames{j}, ...
                'swing');
            end

            % structs for this pData
            thisStanceMatchStepParams.(stepParamNames{j}) = thisStanceStepVal;
            thisSwingMatchStepParams.(stepParamNames{j}) = thisSwingStepVal;
        end

        % get interpolated values for all body parameters
        for j = 1:length(bodytrajParamNames)
            thisMatchBodytrajParams.(bodytrajParamNames{j}) = ...
                interp1(bodytraj.tZeroed, ...
                bodytraj.(bodytrajParamNames{j}), sampTime, 'spline');
        end


        % find all NaNs in PCA variables, note which rows those are found
        % in
        nanRows = [];

        for j = 1:size(pcaVarsMat,1)
            if (any(isnan(pcaVarsMat(j,:))))
                nanRows = [nanRows; j];
            end
        end

        % remove rows containing NaNs from PCA matrix as well as all
        %  interpolated variables
        pcaVarsMat(nanRows,:) = [];

        for j = 1:length(stepParamNames)
            thisStanceMatchStepParams.(stepParamNames{j})(nanRows,:) = [];
            thisSwingMatchStepParams.(stepParamNames{j})(nanRows,:) = [];
        end

        for j = 1:length(bodytrajParamNames)
            thisMatchBodytrajParams.(bodytrajParamNames{j})(nanRows,:) = [];
        end

        
        % add these values to running tracker across pData files
        set1PCVars = [set1PCVars; pcaVarsMat(:,1:set1NumVars)];
        set2PCVars = [set2PCVars; pcaVarsMat(:, ...
            (set1NumVars+1):(set1NumVars + set2NumVars))];

        for j = 1:length(stepParamNames)
            matchStanceStepParams.(stepParamNames{j}) = ...
                [matchStanceStepParams.(stepParamNames{j}); ...
                thisStanceMatchStepParams.(stepParamNames{j})];

            matchSwingStepParams.(stepParamNames{j}) = ...
                [matchSwingStepParams.(stepParamNames{j}); ...
                thisSwingMatchStepParams.(stepParamNames{j})];
        end

        for j = 1:length(bodytrajParamNames)
            matchBodytrajParams.(bodytrajParamNames{j}) = ...
                [matchBodytrajParams.(bodytrajParamNames{j}); ...
                thisMatchBodytrajParams.(bodytrajParamNames{j})];
        end

        % update pData tracker with start and end indices for this pData
        %  file
        pDataFiles.inds = [pDataFiles.ind; ...
            currNumPts + 1, currNumPts + size(pcaVarsMat,1)];
        currNumPts = currNumPts + size(pcaVarsMat,1);

    end

    % remove pData files that didn't contribute any pts
    pDataFiles.names(rmvPDataInd) = [];

    % perform PCA - across all pData files
    [set1.coeff, set1Score, set1.latent, set1.tsquared, ...
        set1.explained, set1.mu] = pca(set1PCVars);

    [set2.coeff, set2Score, set2.latent, set2.tsquared, ...
        set2.explained, set2.mu] = pca(set2PCVars);

    % save output
    saveFileFullPath = [saveFileDir filesep saveFileName '.mat'];

    save(saveFileFullPath, 'set1', 'set2', 'set1PCVars', 'set2PCVars',...
        'set1Score', 'set2Score', 'matchStanceStepParams', ...
        'matchSwingStepParams', 'matchBodytrajParams', ...
        'interpFrameRate', 'legIDs', 'pDataFiles', '-v7.3');

    % print some outputs to screen
    fprintf('Set 1, PC1, variance explained = %.2f%%\n', set1.explained(1));
    fprintf('Set 2, PC1, variance explained = %.2f%%\n', set2.explained(1));
end