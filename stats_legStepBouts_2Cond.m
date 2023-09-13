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
        thisParamVals = squeeze(thisParamVals);

        paramVals{i} = thisParamVals;

    end

    % run t-tests

    allP = zeros(6,1);

%     for i = 1:6
%         [~,allP(i)] = ttest2(paramVals{1}(i,:), paramVals{2}(i,:));
%     end

    % for circular parameter
    for i = 1:6
        paramVal1 = paramVals{1}(i,:);
        paramVal1(isnan(paramVal1)) = [];
        paramVal1 = deg2rad(paramVal1);

        paramVal2 = paramVals{2}(i,:);
        paramVal2(isnan(paramVal2)) = [];
        paramVal2 = deg2rad(paramVal2);

        [allP(i), ~] = circ_wwtest(paramVal1, paramVal2);
    end