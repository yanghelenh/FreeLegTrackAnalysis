% plotScatterPCAStepParamsInterp.m
%
% Function to plot scatterplot of output of savePCAinterpStepParam()
% Plots each step param set 1 PC1 vs set 2 PC1, with option to color by
%  another step parameter or body parameter
% Select output file of savePCAinterpStepParam() through GUI
%
% INPUTS:
%   datDir - path to folder containing savePCAinterpStepParam() output files
%   plotType - string for type of plot to make
%       'binScatter' - binned scatterplot, counts
%       'scatter' - regular scatterplot
%       'stepParam' - scatterplot with points colored by specified matched
%         step parameter value
%       'bodyParam' - scatterplot with points colored by specified matched
%         bodytraj parameter value
%   numBins - for 'binScatter' type plot, number of bins per dimension,
%   whichParam - for stepParam and bodyParam plots, which parameter to plot
%     for bodyParam, string for parameter name. For stepParam, struct with
%     fields for which param (name), which leg (leg), and which phase
%     (phase)
%   xyLine - boolean for whether to plot y = x line
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 7/31/23 - HHY
%
% UPDATED:
%   7/31/23 - HHY
%
function plotScatterPCAStepParamsInterp(datDir, plotType, numBins, ...
    whichParam, xyLine)

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    if (any(strcmpi(plotType, {'binScatter', 'scatter'})))
        load(fullFilePath, 'set1', 'set2', 'set1Score', 'set2Score');
    elseif (strcmpi(plotType, 'bodyParam'))
        load(fullFilePath, 'set1', 'set2', 'set1Score', 'set2Score', ...
            'matchBodytrajParams');
    elseif (strcmpi(plotType, 'stepParam'))
        if (strcmpi(whichParam.phase, 'stance'))
            load(fullFilePath, 'set1', 'set2', 'set1Score', 'set2Score', ...
                'matchStanceStepParams', 'legIDs');
        elseif (strcmpi(whichParam.phase, 'swing'))
            load(fullFilePath, 'set1', 'set2', 'set1Score', 'set2Score', ...
                'matchSwingStepParams', 'legIDs');
        end
    end

    figure;

    % scatterplot all turn bouts in PC space
    if (strcmpi(plotType, 'binScatter'))
        binscatter(set1Score(:,1), set2Score(:,1), numBins);
    elseif (strcmpi(plotType, 'bodyParam'))
        thisBodyParam = matchBodytrajParams.(whichParam);
        scatter(set1Score(:,1), set2Score(:,1), [], ...
            thisBodyParam, 'filled', ...
            'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0.03);
        colormap(gca, 'cool');
        c = colorbar;
%         caxis([-300 300]);
        caxis([-0.01 0.04]);
        c.Label.String = whichParam;

    elseif (strcmpi(plotType, 'stepParam'))
        if (strcmpi(whichParam.phase, 'stance'))
            thisStepParam = matchStanceStepParams.(whichParam.name);
        elseif (strcmpi(whichParam.phase, 'swing'))
            thisStepParam = matchSwingStepParams.(whichParam.name);
        end

        thisLegInd = legIDs.ind(strcmpi(whichParam.leg, legIDs.name));
        thisVal = thisStepParam(:,thisLegInd);

        scatter(set1.score(:,1), set2.score(:,1), [], ...
            thisVal, 'filled', ...
            'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
        colormap(gca, 'cool');
        c = colorbar;
        c.Label.String = whichParam.name;
    else
        scatter(set1Score(:,1), set2Score(:,1), 'filled', ...
            'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0.03);
    end

    hold on;

    % plot reference lines
    xLimits = xlim;
    yLimits = ylim;

    minPt = min([xLimits(1) yLimits(1)]);
    maxPt = max([xLimits(2) yLimits(2)]);

    if (xyLine)
        plot([minPt maxPt], [minPt maxPt], 'k', 'LineWidth', 2);
    end

    % plot x = 0 line
    line([0 0], [minPt maxPt], 'Color', 'k');
    % plot y = 0 line
    line([minPt maxPt], [0 0], 'Color', 'k');

    xAxisLabel = sprintf('%s, PC1: %.2f', set1.name, set1.explained(1));
    xlabel(xAxisLabel);

    yAxisLabel = sprintf('%s, PC1: %.2f', set2.name, set2.explained(1));
    ylabel(yAxisLabel);

    xlim(xLimits);
    ylim(yLimits);
end