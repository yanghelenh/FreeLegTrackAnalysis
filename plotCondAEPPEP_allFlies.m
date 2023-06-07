% plotCondAEPPEP_allFlies.m
%
% Function that takes output files from saveLegStepParamCond_indpt() and
%  generates 2 plots: AEP and PEP
% Either plots mean +/- std dev for each fly or across all flies.
% Can select more than 1 output file, will plot on same graph
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   xyScale - scale for plots, as [min max]
%   indivFlies - boolean for whether to plot individual flies
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 5/9/23 - HHY
%
% UPDATED:
%   5/9/23 - HHY
%
function plotCondAEPPEP_allFlies(datDir, xyScale, indivFlies)

    numLegs = 6;

    % prompt user to select output files from saveLegStepParamCond_indpt()
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

    % initialize figures
    aepFig = figure;
    pepFig = figure;
    c = colormap('lines');

    % preallocate legend str
    legendStr = cell(numFiles,1);

    % preallocate
    AEPxMeans = zeros(numLegs, numFiles);
    AEPyMeans = zeros(numLegs, numFiles);
    PEPxMeans = zeros(numLegs, numFiles);
    PEPyMeans = zeros(numLegs, numFiles);
    
    AEPxSDs = zeros(numLegs, numFiles);
    AEPySDs = zeros(numLegs, numFiles);
    PEPxSDs = zeros(numLegs, numFiles);
    PEPySDs = zeros(numLegs, numFiles);
    
    ns = zeros(numFiles, numLegs);

    for i = 1:numFiles

        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];


        % load data - stance only for AEP/PEP
        load(outputFullPath, 'selStanceParams', 'legIDs', 'cond');

        % check if individual flies are being plotted or not
        if (indivFlies)
            % WRITE THIS LATER %
            % 5/9/23 %
        else
%             % compute mean for each leg
%             numLegs = length(legIDs.ind);


            % loop through each leg
            for j = 1:numLegs
                thisLeg = legIDs.ind(j);
                thisLegLog = selStanceParams.stepWhichLeg == thisLeg;

                AEPxMeans(j,i) = mean(selStanceParams.stepAEPX(thisLegLog));
                AEPyMeans(j,i) = mean(selStanceParams.stepAEPY(thisLegLog));
                PEPxMeans(j,i) = mean(selStanceParams.stepPEPX(thisLegLog));
                PEPyMeans(j,i) = mean(selStanceParams.stepPEPY(thisLegLog));
    
                AEPxSDs(j,i) = std(selStanceParams.stepAEPX(thisLegLog));
                AEPySDs(j,i) = std(selStanceParams.stepAEPY(thisLegLog));
                PEPxSDs(j,i) = std(selStanceParams.stepPEPX(thisLegLog));
                PEPySDs(j,i) = std(selStanceParams.stepPEPY(thisLegLog));

                ns(i,j) = sum(thisLegLog);
            end



            % get legend str, named based on condition
            thisLegendStr = [];
            for j = 1:length(cond.whichParam)
                thisLegendStr = [thisLegendStr cond.whichParam{j}...
                    cond.cond{j}];
                if (j~=length(cond.whichParam))
                    thisLegendStr = [thisLegendStr ', '];
                end
            end

            % add n to legendStr
            thisLegendStr = [thisLegendStr '; n = ' num2str(ns(i,:))];

            legendStr{i} = thisLegendStr;

        end
    end

    % plot this file's data into AEP figure
    figure(aepFig);

    for i = 1:numFiles
        errorbar(AEPyMeans(:,i), AEPxMeans(:,i), AEPxSDs(:,i), ...
            AEPxSDs(:,i), AEPySDs(:,i), AEPySDs(:,i), ...
            'Marker','x', 'LineStyle','none','Color',c(i,:), ...
            'LineWidth',1.5);
    
        hold on;
    end

    title('AEP');
    % x and y lims
    xlim(xyScale);
    ylim(xyScale);

    axis('equal');

    % x and y labels
    xlabel('Body Lengths <-L - R->');
    ylabel('Body Lengths <-P - A->')

    % reverse y axis (x values) so head (neg vals) is at top
    set(gca, 'YDir','reverse');

    legend(legendStr); 

    

    % plot this file's data into PEP figure
    figure(pepFig);

    for i = 1:numFiles
        errorbar(PEPyMeans(:,i), PEPxMeans(:,i), PEPxSDs(:,i), ...
            PEPxSDs(:,i), PEPySDs(:,i), PEPySDs(:,i), ...
            'Marker','x', 'LineStyle','none','Color',c(i,:), ...
            'LineWidth',1.5);
    
        hold on;
    end

    title('PEP');

    % x and y lims
    xlim(xyScale);
    ylim(xyScale);

    axis('equal');

    % x and y labels
    xlabel('Body Lengths <-L - R->');
    ylabel('Body Lengths')

    % reverse y axis (x values) so head (neg vals) is at top
    set(gca, 'YDir','reverse');

    legend(legendStr);
end

