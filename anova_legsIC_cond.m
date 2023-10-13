% anova_legsIC_cond.m
%
% Function to compute 2-way ANOVA, leg identity as
%  one factor and the condition identity as the other factor (select 2 or
%  more output files, each represents one condition).
% Also performs multiple comparisons
% Runs on output files of saveLegStepParamCond_bouts()
%
% INPUTS:
%
% OUTPUTS:
%
% CREATED: 9/12/23 - HHY
%
% UPDATED:
%   9/12/23 - HHY
%
function anova_legsIC_cond(datDir, whichParam, whichPhase, saveName, saveDir)

    circStepParams = {'stepDirections'};
    
    disp('Select output files from saveLegStepParamCond_bouts()');
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select cond_bouts files', datDir, 'MultiSelect', 'on');

    if (~iscell(outputFNames))
        disp('Function requires at least 2 conditions');
        return;
    end

    numCond = length(outputFNames);

    % prompt user to select reference file from 
    %  saveLegStepParamCond_indpt()
    disp('Select reference file');
    [refOutputName, refOutputPath] = uigetfile('*.mat', ...
        'Select reference file', datDir, 'MultiSelect', 'off');

    % load reference
    refOutputFullPath = [refOutputPath filesep refOutputName];

    if strcmpi(whichPhase, 'stance') % stance
        load(refOutputFullPath, 'selStanceParams', 'legIDs');
        refParamVal = selStanceParams.(whichParam);
    else % swing
        load(refOutputFullPath, 'selSwingParams', 'legIDs');
        refParamVal = selSwingParams.(whichParam);
    end

    % number of legs
    numLegs = length(legIDs.ind);

    % preallocate
    refMeans = zeros(numLegs, 1);

    % loop through each leg
    for j = 1:numLegs
        thisLeg = legIDs.ind(j);
        if strcmpi(whichPhase, 'stance') % stance
            thisLegLog = selStanceParams.stepWhichLeg == thisLeg;
        else
            thisLegLog = selSwingParams.stepWhichLeg == thisLeg;
        end
        
        % compute mean, std; check if circular param
        if(any(strcmpi(whichParam, circStepParams)))
            % convert to radians
            refParamValRad = deg2rad(wrapTo180(rmoutliers(wrapTo360(refParamVal(thisLegLog)))));
            % get mean
            refMeans(j) = rad2deg(circ_mean(refParamValRad));
        else
            refMeans(j) = mean(rmoutliers(refParamVal(thisLegLog)));
        end
    end

    % initialize - vectors to make table
    legIC = [];
    cond = [];
    neuromere = [];
    paramVal = [];

    for i = 1:numCond
        outputFullPath = [outputPath filesep outputFNames{i}];

        condName = string(outputFNames{i}(1:(end-4)));

        if strcmpi(whichPhase, 'stance')
            load(outputFullPath, 'selStanceParams');
            thisParamVals = selStanceParams.(whichParam);

        else
            load(outputFullPath, 'selSwingParams');
            thisParamVals = selSwingParams.(whichParam);
        end

        % reshape as column vector, with first half as ipsi, second as contra
        thisParamVals = squeeze(thisParamVals);
        for j = 1:numLegs
            thisParamVals(j,:) = thisParamVals(j,:) - refMeans(j);
        end

        thisParamVals = thisParamVals';
        thisParamVals = thisParamVals(:);

        % wrap to 360 for stance step directions 
        if any(strcmpi(whichParam,circStepParams)) && strcmpi(whichPhase,'stance')
            thisParamVals = wrapTo360(thisParamVals);
        end

        % append to output
        paramVal = [paramVal; thisParamVals];

        % number of samples
        numSamps = length(thisParamVals);
        
        % get vector of ipsi/contra 
%         thisLegIC = [repmat("Ipsi",numSamps/2,1); ...
%             repmat("Contra", numSamps/2, 1)];

%         thisNeuromere = [repmat("Front",numSamps/6,1); repmat("Mid", numSamps/6, 1);
%             repmat("Hind", numSamps/6, 1); repmat("Front", numSamps/6, 1);
%             repmat("Mid", numSamps/6, 1); repmat("Hind", numSamps/6, 1);];

        thisLegIC = [repmat("R1",numSamps/6,1); repmat("R2", numSamps/6, 1);
            repmat("R3", numSamps/6, 1); repmat("L1", numSamps/6, 1);
            repmat("L2", numSamps/6, 1); repmat("L3", numSamps/6, 1);];
        % append
        legIC = [legIC; thisLegIC];
        neuromere = [neuromere; thisNeuromere];

        % get vector of cond
        thisCond = repmat(condName,numSamps,1);
        cond = [cond; thisCond];
    end


    % perform unbalanced ANOVA
%     [p, tbl, stats] = anovan(paramVal, {legIC, neuromere, cond}, ...
%         'model','interaction');
    [p, tbl, stats] = anovan(paramVal, {legIC, cond}, ...
        'model','interaction');

    % multiple comparison tests
    [multCmpResults,~,~,gNames] = multcompare(stats, 'Dimension',[1 2], ...
        'CType', 'tukey-kramer');

    % save output
    saveFullPath = [saveDir filesep saveName '.mat'];

    save(saveFullPath, 'whichParam', 'whichPhase', 'p', 'tbl', 'stats', ...
        'multCmpResults', 'gNames');

end