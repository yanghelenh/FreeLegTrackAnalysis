% saveContLegStepParamCond_indpt.m
%
% Function that saves continuous leg step parameters for all steps that
%  meet the user specified conditions on bodytraj
% Treats all time points as independent, as long as bodytraj condition(s) 
%  are met, unlike saveContLegStepParamCond_bouts.m
% Modification of saveLegStepParamCond_indpt(), which operates on steps,
%  while this function operates on continuous leg parameters
% User selects one or more pData files through GUI or input
% 
% INPUTS:
%   cond - struct of conditions, if multiple conditions, treats it as AND
%       whichParam - cell array (even if 1 element) on which bodytraj field
%           to condition on, one for each condition
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       selLegStepsCont - all time points, all pData files, that meet 
%           conditions. Each field is different step parameter. Also
%           includes field that indicates which leg
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to indices in legSteps) that
%               belong to each pData file, as start and end indices
%       legIDs - legStepsCont.legIDs
%       cond - same as INPUT
%
% CREATED: 9/6/23 - HHY
%
% UPDATED:
%   9/6/23 - HHY
%
function saveContLegStepParamCond_indpt(cond, pDataFNames, pDataPath, ...
    saveFilePath, saveFileName)

    % names of all step parameters to save
    stepParamNames = {'AEPX', 'PEPX', 'AEPY', 'PEPY', 'stepLengthX', ...
        'stepLengthY', 'stepLength', 'stepDirection'};

    % all the step parameters where values need to be * -1 for left turns
    flipStepParams = {'AEPY', 'PEPY', 'stepLengthY', 'stepDirection'};

    % bodytraj parameters that are L/R asymmetric contain one of the
    %  following
    lrAsymBTParams = {'angVel', 'latVel'};

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

    % preallocate 
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];

    for i = 1:length(stepParamNames)
        selLegStepsCont.(stepParamNames{i}) = [];
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

        % check if this pData file has legStepsCont, bodytraj, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'bodytraj')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            continue;
        end

        % load legTrack and moveNotMove from pData
        load(pDataFullPath, 'moveNotMove', 'bodytraj', ...
            'legStepsCont');

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

        % check if we're conditioning on any L/R asymmetric parameters
        % if yes, allocate left and right turns separately
        if(any(contains(cond.whichParam, lrAsymBTParams))) 
            sepLR = true; % logical for whether we're separating L/R
            % preallocate logical for valid indices of bodytraj
            btRightValInd = true(size(bodytraj.tZeroed));
            btLeftValInd = true(size(bodytraj.tZeroed));

            % initialize - keep track of indices of steps to include
            inclStepIndRight = [];
            inclStepIndLeft = [];
        else % if no, no need to separate left and right
            sepLR = false;
            % preallocate logical for valid indices of bodytraj
            btValInd = true(size(bodytraj.tZeroed));

            % initialize - keep track of indices of steps to include
            inclStepInd = [];
        end

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

            % if we're separating left and right (at least 1 param is
            %  asymmetric)
            if(sepLR)
                % get valid indices for right turns
                thisValIndRight = eval(['thisCondVar' cond.cond{j}]);
                % combine with all other conditions
                btRightValInd = btRightValInd & thisValIndRight;

                % if this is a parameter to be flipped, flip it for left
                if(contains(cond.whichParam{j},lrAsymBTParams))
                    thisCondVar = thisCondVar * -1;
                end

                % left turns
                thisValIndLeft = eval(['thisCondVar' cond.cond{j}]);
                btLeftValInd = btLeftValInd & thisValIndLeft;
            else % not separating left and right
                % bodytraj indices where this condition is true
                thisValInd = eval(['thisCondVar' cond.cond{j}]);
    
                % combine with all other conditions
                btValInd = btValInd & thisValInd;
            end
        end

        % filter for moving bouts only
        if(sepLR)
            btRightValInd(moveNotMove.notMoveInd) = false;
            btLeftValInd(moveNotMove.notMoveInd) = false;
        else
            btValInd(moveNotMove.notMoveInd) = false;
        end



        % get start index
        trialStartInd = size(selLegStepsCont.(stepParamNames{1}), 1) + 1;

        % loop through all step parameters, use valid indices to get output
        for j = 1:length(stepParamNames)
            % if separating left and right
            if (sepLR)
                if (~isempty(btRightValInd) || ~isempty(btLeftValInd))
                    % right
                    thisRightParamVal = ...
                        legStepsCont.(stepParamNames{j})(btRightValInd,:);
                    % left - swap left and right
                    thisLeftParamVal = ...
                        legStepsCont.(stepParamNames{j})(btLeftValInd, ...
                        matchedLegInd);
    
                    % if this is a parameter where the value needs to be
                    % inverted
                    if any(strcmpi(stepParamNames{j}, flipStepParams))
                        thisLeftParamVal = thisLeftParamVal * -1;
                    end
    
                    % combine left and right
                    thisParamVal = [thisRightParamVal; thisLeftParamVal];

                end
                
            % not separating left and right
            else
                if (~isempty(btValInd))
                    thisParamVal = ...
                        legStepsCont.(stepParamNames{j})(btValInd,:);
                end
            end

            % concatenate
            if ~isempty(thisParamVal)
                selLegStepsCont.(stepParamNames{j}) = [...
                    selLegStepsCont.(stepParamNames{j}); thisParamVal];
            end
        end

        % end index for steps for this trial
        trialEndInd = size(selLegStepsCont.(stepParamNames{1}), 1);

        % check if this trial contributes any time points
        if (trialEndInd >= trialStartInd)
            thisTrialInds = [trialStartInd trialEndInd];
    
            pDataFiles.inds = [pDataFiles.inds; thisTrialInds];
        
        % otherwise, remove from pData list
        else
            rmvInd = [rmvInd; i];
        end
    end

    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    % save output file
    fullSavePath = [saveFilePath filesep saveFileName '.mat'];

    save(fullSavePath, 'selLegStepsCont', 'pDataFiles', 'cond', ...
        'legIDs', '-v7.3');
end