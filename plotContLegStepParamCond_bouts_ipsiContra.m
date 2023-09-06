% plotContLegStepParamCond_bouts_ipsiContra.m
%
% Function to plot overlay of ipsi and contra legs for output of 
%  saveContLegStepParamCond_bouts(), 1 file.
%
% INPUTS:
%   whichParam - name of parameter to plot
%   datDir - directory with output files
%   yScale - y scale for plots, as [min max]
%   semError - boolean for whether to plot SEM; if false, plots std dev
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/6/23 - HHY
%
% UPDATED:
%   9/6/23 - HHY
%
function plotContLegStepParamCond_bouts_ipsiContra(whichParam, datDir, ...
    yScale, semError)


    % all step parameters that are circular variables - need to use
    %  circular stats
    circStepParams = {'stepDirections'};

    % prompt user to select output files from saveContLegStepParamCond_bouts()
    disp('Select output file from saveContLegStepParamCond_bouts()');
    [outputFName, outputPath] = uigetfile('*.mat', ...
        'Select cond_bouts files', datDir, 'MultiSelect', 'off');
        
    outputFullPath = [outputPath outputFName];

    % load data - with appropriate error
    if (semError) % SEM
        load(outputFullPath, 'selLegStepsCont', 'legStepsContMean', ...
            'legStepsContSEM', 'legT', 'cond', 'fwdVelCond');
        paramErr = legStepsContSEM.(whichParam);
    else
        load(outputFullPath, 'selLegStepsCont', 'legStepsContMean', ...
            'legStepsContStd', 'legT', 'cond', 'fwdVelCond');
        paramErr = legStepsContStd.(whichParam);
    end


    % if circular parameter, wrap to 360 for plotting
    if(any(strcmpi(whichParam, circStepParams)))
        paramMean = wrapTo360(legStepsContMean.(whichParam));
    else % not circular
        paramMean = legStepsContMean.(whichParam);
    end

    numBouts = size(selLegStepsCont.(whichParam), 1);


    % initialize figure
    f = figure;
    c = colormap('lines');

    for i = 1:3
        subplot(3,1,i);
        hold on;

        % plot
        errorbar(legT', paramMean(:,i), paramErr(:,i), ...
            'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(1,:));

        errorbar(legT', paramMean(:,i+3), paramErr(:,i+3), ...
            'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(2,:));

        % axis scale and label

        % for AEP
%         if (i == 1)
%             ylim([-0.9 -0.5]);
%         elseif (i == 2)
%             ylim([-0.3 0.1]);
%         else
%             ylim([0.15 0.55]);
%         end

        % for PEP
%         if (i == 1)
%             ylim([-0.6 -0.2]);
%         elseif (i == 2)
%             ylim([0 0.4]);
%         else
%             ylim([0.45 0.85]);
%         end

        % line at t = 0
        line([0 0], yScale, 'LineWidth', 1, 'Color', 'k');

        ylim(yScale);

        xScale = [legT(1), legT(end)];
        xlim(xScale);

        xlabel('Time (s)');
        ylabel(whichParam);

        % label legs
        if (i == 1)
            title('Front Legs');
        elseif (i == 2)
            title('Mid legs');
        else
            title('Hind legs');
        end

        % legend
        legend({'Ipsi', 'Contra'});
    end

    sgtitle(whichParam);

end