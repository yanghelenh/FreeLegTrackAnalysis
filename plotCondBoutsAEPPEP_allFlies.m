% plotCondBoutsAEPPEP_allFlies.m
%
% Function that takes output files from saveLegStepParamCond_bouts() and
%  generates 2 plots (AEP and PEP) for each time step 
% Plots as either mean +/- SEM or mean +/- std dev
% Can select more than 1 output file, will plot on same graph
% Has option to plot reference steps (from saveLegStepParamCond_indpt()).
%  Select one file through GUI
% Select output files through GUI
% Note: make sure all output files selected have same maxNumSteps. Plotting
%  will get confused otherwise
%
% INPUTS:
%   datDir - directory with output files
%   xyScale - scale for plots, as [min max]
%   semError - boolean for whether to plot SEM; if false, plots std dev
%   ref - boolean for whether to plot reference
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 5/15/23 - HHY
%
% UPDATED:
%   5/15/23 - HHY
%   6/5/23 - HHY - adds the option to include a reference file
%   6/7/23 - HHY - removes outliers from reference file before computing
%       mean and std
%
function plotCondBoutsAEPPEP_allFlies(datDir, xyScale, semError, ref)

    % prompt user to select output files from saveLegStepParamCond_bouts()
    disp('Select output files from saveLegStepParamCond_bouts()');
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select AEP PEP files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFiles = length(outputFNames);
    else
        numFiles = 1;
    end

    % preallocate
    allAEPXMeans = cell(numFiles, 1);
    allAEPYMeans = cell(numFiles, 1);
    allPEPXMeans = cell(numFiles, 1);
    allPEPYMeans = cell(numFiles, 1);
    allAEPXErrs = cell(numFiles, 1);
    allAEPYErrs = cell(numFiles, 1);
    allPEPXErrs = cell(numFiles, 1);
    allPEPYErrs = cell(numFiles, 1);
    numBoutsAEP = cell(numFiles, 1);
    numBoutsPEP = cell(numFiles, 1);
    legendStr = cell(numFiles, 1);

    % if reference
    if (ref)
        % prompt user to select reference file from 
        %  saveLegStepParamCond_indpt()
        disp('Select reference file');
        [refOutputName, refOutputPath] = uigetfile('*.mat', ...
            'Select reference AEP PEP file', datDir, 'MultiSelect', 'off');

        % load reference
        refOutputFullPath = [refOutputPath filesep refOutputName];
        load(refOutputFullPath, 'selStanceParams', 'legIDs');

        % number of legs
        numLegs = length(legIDs.ind);

        % preallocate
        refAEPxMeans = zeros(numLegs, 1);
        refAEPyMeans = zeros(numLegs, 1);
        refPEPxMeans = zeros(numLegs, 1);
        refPEPyMeans = zeros(numLegs, 1);
        
        refAEPxErrs = zeros(numLegs, 1);
        refAEPyErrs = zeros(numLegs, 1);
        refPEPxErrs = zeros(numLegs, 1);
        refPEPyErrs = zeros(numLegs, 1);

        refN = zeros(1, numLegs);
    end

    % loop through all files, save AEP and PEP matrices for each file
    % save also, cond, numBouts
    for i = 1:numFiles
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data - stance only for AEP/PEP
        % load appropriate error variable
        if (semError)
            load(outputFullPath, 'stanceParamMeans', 'stanceParamSEM', ...
                'stanceParamN', 'maxNumSteps', 'cond');

            allAEPXErrs{i} = stanceParamSEM.stepAEPX;
            allAEPYErrs{i} = stanceParamSEM.stepAEPY;
            allPEPXErrs{i} = stanceParamSEM.stepPEPX;
            allPEPYErrs{i} = stanceParamSEM.stepPEPY;
        else
            load(outputFullPath, 'stanceParamMeans', 'stanceParamStd', ...
                'stanceParamN', 'maxNumSteps', 'cond');

            allAEPXErrs{i} = stanceParamStd.stepAEPX;
            allAEPYErrs{i} = stanceParamStd.stepAEPY;
            allPEPXErrs{i} = stanceParamStd.stepPEPX;
            allPEPYErrs{i} = stanceParamStd.stepPEPY;
        end

        % other variables
        allAEPXMeans{i} = stanceParamMeans.stepAEPX;
        allAEPYMeans{i} = stanceParamMeans.stepAEPY;
        allPEPXMeans{i} = stanceParamMeans.stepPEPX;
        allPEPYMeans{i} = stanceParamMeans.stepPEPY;
        numBoutsAEP{i} = stanceParamN.stepAEPX;
        numBoutsPEP{i} = stanceParamN.stepPEPX;

        % get legend str, named based on condition
        thisLegendStr = [];
        for j = 1:length(cond.whichParam)
            thisLegendStr = [thisLegendStr cond.whichParam{j}...
                cond.cond{j}];
            if (j~=length(cond.whichParam))
                thisLegendStr = [thisLegendStr ', '];
            end
        end
        legendStr{i} = thisLegendStr;
    end

    % if reference, compute means and errors
    if (ref)
        % loop through each leg
        for j = 1:numLegs
            thisLeg = legIDs.ind(j);
            thisLegLog = selStanceParams.stepWhichLeg == thisLeg;
    
            % compute means, remove outliers
            refAEPxMeans(j) = mean(rmoutliers(selStanceParams.stepAEPX(thisLegLog)));
            refAEPyMeans(j) = mean(rmoutliers(selStanceParams.stepAEPY(thisLegLog)));
            refPEPxMeans(j) = mean(rmoutliers(selStanceParams.stepPEPX(thisLegLog)));
            refPEPyMeans(j) = mean(rmoutliers(selStanceParams.stepPEPY(thisLegLog)));
    
            % n for this leg
            thisN = length(rmoutliers(selStanceParams.stepAEPX(thisLegLog)));
            refN(j) = thisN;

            % compute errors, remove outliers
            if (semError) % SEM
                refAEPxErrs(j) = std(rmoutliers(selStanceParams.stepAEPX(thisLegLog))) ...
                    / sqrt(thisN);
                refAEPyErrs(j) = std(rmoutliers(selStanceParams.stepAEPY(thisLegLog)))...
                    / sqrt(thisN);
                refPEPxErrs(j) = std(rmoutliers(selStanceParams.stepPEPX(thisLegLog)))...
                    / sqrt(thisN);
                refPEPyErrs(j) = std(rmoutliers(selStanceParams.stepPEPY(thisLegLog)))...
                    / sqrt(thisN);
            else % std dev
                refAEPxErrs(j) = std(rmoutliers(selStanceParams.stepAEPX(thisLegLog)));
                refAEPyErrs(j) = std(rmoutliers(selStanceParams.stepAEPY(thisLegLog)));
                refPEPxErrs(j) = std(rmoutliers(selStanceParams.stepPEPX(thisLegLog)));
                refPEPyErrs(j) = std(rmoutliers(selStanceParams.stepPEPY(thisLegLog)));
            end
    
        end
    end

    % get step time points (for x-axis): 0 for step at peak
    stepTPts = [fliplr((1:maxNumSteps) * -1), 0, 1:maxNumSteps];
    numStepTPts = length(stepTPts);

    % loop through all steps
    for i = 1:numStepTPts

        % initialize figures
        aepFig = figure;
        pepFig = figure;
        c = colormap('lines');
    
        % preallocate legend str
        if(ref)
            legendStrAEP = cell(numFiles + 1,1);
            legendStrPEP = cell(numFiles + 1,1);
        else
            legendStrAEP = cell(numFiles,1);
            legendStrPEP = cell(numFiles,1);
        end

        % AEP plot
        figure(aepFig);

        % if reference, plot those AEPs in black
        if(ref)
            errorbar(refAEPyMeans, refAEPxMeans, ...
                refAEPxErrs, refAEPxErrs, ...
                refAEPyErrs, refAEPyErrs, ...
                'Marker','x', 'LineStyle','none','Color','black', ...
                'LineWidth',1.5);
        
            % legend
            thisLegendStr = ['Reference; n = ' num2str(refN)];
    
            legendStrAEP{1} = thisLegendStr;
    
            hold on;
        end

        % loop through all files
        for j = 1:numFiles
            % plot
            errorbar(allAEPYMeans{j}(i,:), allAEPXMeans{j}(i,:), ...
                allAEPXErrs{j}(i,:), allAEPXErrs{j}(i,:), ...
                allAEPYErrs{j}(i,:), allAEPYErrs{j}(i,:), ...
                'Marker','x', 'LineStyle','none','Color',c(j,:), ...
                'LineWidth',1.5);
    
            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j} '; n = ' ...
                num2str(numBoutsAEP{j}(i,:))];

            if(ref)
                legendStrAEP{j + 1} = thisLegendStr;
            else
                legendStrAEP{j} = thisLegendStr;
            end
        end

        % get title string
        ttlStr = ['AEP, Step Num = ' num2str(stepTPts(i))];
        title(ttlStr);

        % x and y lims
        xlim(xyScale);
        ylim(xyScale);
    
        axis('equal');
    
        % x and y labels
        xlabel('Body Lengths <-L - R->');
        ylabel('Body Lengths <-P - A->')
    
        % reverse y axis (x values) so head (neg vals) is at top
        set(gca, 'YDir','reverse');
    
        legend(legendStrAEP); 


        % PEP plot
        figure(pepFig);

        % if reference, plot those AEPs in black
        if(ref)
            errorbar(refPEPyMeans, refPEPxMeans, ...
                refPEPxErrs, refPEPxErrs, ...
                refPEPyErrs, refPEPyErrs, ...
                'Marker','x', 'LineStyle','none','Color','black', ...
                'LineWidth',1.5);
    
            % legend
            thisLegendStr = ['Reference; n = ' num2str(refN)];
    
            legendStrPEP{1} = thisLegendStr;
    
            hold on;
        end

        % loop through all files
        for j = 1:numFiles
            % plot
            errorbar(allPEPYMeans{j}(i,:), allPEPXMeans{j}(i,:), ...
                allPEPXErrs{j}(i,:), allPEPXErrs{j}(i,:), ...
                allPEPYErrs{j}(i,:), allPEPYErrs{j}(i,:), ...
                'Marker','x', 'LineStyle','none','Color',c(j,:), ...
                'LineWidth',1.5);
    
            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j} '; n = ' ...
                num2str(numBoutsPEP{j}(i,:))];

            if(ref)
                legendStrPEP{j+1} = thisLegendStr;
            else
                legendStrPEP{j} = thisLegendStr;
            end
        end

        % get title string
        ttlStr = ['PEP, Step Num = ' num2str(stepTPts(i))];
        title(ttlStr);

        % x and y lims
        xlim(xyScale);
        ylim(xyScale);
    
        axis('equal');
    
        % x and y labels
        xlabel('Body Lengths <-L - R->');
        ylabel('Body Lengths <-P - A->')
    
        % reverse y axis (x values) so head (neg vals) is at top
        set(gca, 'YDir','reverse');
    
        legend(legendStrPEP); 
    end

end