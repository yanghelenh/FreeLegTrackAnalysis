% quick script mod of anova_legsIC_cond , for circular parameters   

circStepParams = {'stepDirections'};

    legNames = {"R1", "R2", "R3", "L1", "L2", "L3"};
    
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
    refVals = {};

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
            refVals{j} = wrapTo180(rmoutliers(wrapTo360(refParamVal(thisLegLog))));
        else
            refMeans(j) = mean(rmoutliers(refParamVal(thisLegLog)));
            refVals{j} = rmoutliers(refParamVal(thisLegLog));
        end
    end

    % initialize - vectors to make table
    legIC = [];
    cond = [];
%     neuromere = [];
    paramVal = [];

    if ~calcDiff
        % add reference
        for i = 1:numLegs
            thisLegIC = repmat(legNames{i},length(refVals{i}),1);
            legIC = [legIC; thisLegIC];
            thisCond = repmat("Ref",length(refVals{i}),1);
            cond = [cond; thisCond];
            paramVal = [paramVal; refVals{i}];
        end
    end

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
            if ~(any(strcmpi(whichParam, circStepParams)))
                if calcDiff
                    thisParamVals(j,:) = thisParamVals(j,:) - refMeans(j);
                else
                    thisParamVals(j,:) = thisParamVals(j,:);
                end
            else
                if calcDiff
                    thisParamVals(j,:) = deg2rad(wrapTo180(thisParamVals(j,:))) - refMeans(j);
                else
                    thisParamVals(j,:) = deg2rad(wrapTo180(thisParamVals(j,:)));
                end
            end

        end

        thisParamVals = thisParamVals';
        thisParamVals = thisParamVals(:);

%         % wrap to 360 for stance step directions 
%         if any(strcmpi(whichParam,circStepParams)) && strcmpi(whichPhase,'stance')
%             thisParamVals = wrapTo360(thisParamVals);
%         end

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
%         neuromere = [neuromere; thisNeuromere];

        % get vector of cond
        thisCond = repmat(condName,numSamps,1);
        cond = [cond; thisCond];
    end


    % perform unbalanced ANOVA
%     [p, tbl, stats] = anovan(paramVal, {legIC, neuromere, cond}, ...
%         'model','interaction');
    if ~(any(strcmpi(whichParam, circStepParams)))
        [p, tbl, stats] = anovan(paramVal, {legIC, cond}, ...
            'model','interaction');
    
        % multiple comparison tests
        [multCmpResults,~,~,gNames] = multcompare(stats, 'Dimension',[1 2], ...
            'CType', 'tukey-kramer');
    else
        nanLog = isnan(paramVal);
        paramVal(nanLog) = [];
        legIC(nanLog) = [];
        cond(nanLog) = [];

        [p, stats] = circ_hktest(paramVal, legIC, cond, 1, {'Leg','Cond'});

    end

    % post-hoc tests
    allCondNames = unique(cond);
    allCondNames(1) = [];

    whichLeg = 3;
    whichCond = 1;

    thisCondLog = strcmpi(legIC,legNames{whichLeg}) & strcmpi(cond,allCondNames(whichCond));
    condVal = paramVal(thisCondLog);

    thisRefLog = strcmpi(legIC,legNames{whichLeg}) & strcmpi(cond,"Ref");
    refVal = paramVal(thisRefLog);

    [p, stats] = circ_wwtest(refVal, condVal);

