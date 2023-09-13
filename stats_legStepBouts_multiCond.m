whichParam = 'stepDirections';
whichPhase = 'stance';


    disp('Select output files from saveLegStepParamCond_bouts()');
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select cond_bouts files', datDir, 'MultiSelect', 'on');

    numCond = length(outputFNames);

    paramVals = {};

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
%         thisParamVals = squeeze(thisParamVals);
        thisParamVals = squeeze(thisParamVals);

        paramVals{i} = thisParamVals;

    end

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
        refParamLegs = selStanceParams.stepWhichLeg;
    else % swing
        load(refOutputFullPath, 'selSwingParams', 'legIDs');
        refParamVal = selSwingParams.(whichParam);
        refParamLegs = selSwingParams.stepWhichLeg;
    end

    % run anova

    allP = zeros(6,1);
    allStats = {};
    allC = {};

    for i = 1:6
        refVal = refParamVal(refParamLegs == i);
        allParamVals = refVal;
        allParamCond = zeros(length(refVal),1);
        for j = 1:numCond
            thisParamVal = paramVals{j}(i,:);
            thisParamVal(isnan(thisParamVal)) = [];
            thisParamVal = deg2rad(thisParamVal);
            thisParamVal = thisParamVal';
%             allParamVals = [allParamVals; paramVals{j}(i,:)'];
            allParamVals = [allParamVals; thisParamVal];
%             allParamCond = [allParamCond; repmat(j,size(paramVals{j},2),1)];
            allParamCond = [allParamCond; repmat(j,length(thisParamVal),1)];
        end

        [allP(i), allStats{i}] = circ_wwtest(deg2rad(allParamVals),allParamCond);
%         [allP(i),~,allStats{i}] = anova1(allParamVals,allParamCond);
% 
%         allC{i} = multcompare(allStats{i});
    end


