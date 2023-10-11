% plotNormContLegStepParamCond_bouts_ipsiContra.m
%
% Normalize each bout's step param by subtracting out min and dividing by
%  range. Compute FWHM
%
% INPUTS:
%   whichParam - name of parameter to plot
%   datDir - directory with output files
%   yScale - y scale for plots, as [min max]
%   semError - boolean for whether to plot SEM; if false, plots std dev
%
% OUTPUTS:
%   p - p-values for paired t-test for FWHM 

function p = plotNormContLegStepParamCond_bouts_ipsiContra(whichParam, datDir,...
    fwhmValInd, yScale, semError)

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

    % get this param val
    paramVal = selLegStepsCont.(whichParam);

    % preallocate
    paramValNorm = zeros(size(paramVal));
    fwhmAll = zeros(size(paramVal,2),size(paramVal,3));

    % loop through all legs
    for i = 1:6
        thisLegParamVal = squeeze(paramVal(:,i,:));
        % loop through all bouts
        for j = 1:size(thisLegParamVal,2)
            thisParamVal = thisLegParamVal(:,j);
            % normalize to 0-1 by subtracting min and dividing by range
            thisNormParamVal = (thisParamVal - min(thisParamVal)) / ...
                (max(thisParamVal) - min(thisParamVal));

            % save to output matrix
            paramValNorm(:,i,j) = thisNormParamVal;

            % get fwhm
            if (i < 4)
                thisFwhm = fwhm(legT, thisNormParamVal,fwhmValInd,'min');
            else
                thisFwhm = fwhm(legT, thisNormParamVal,fwhmValInd,'max');
            end

            fwhmAll(i,j) = thisFwhm;
        end
    end

    % get mean, std, sem of normalized param val
    % preallocate
    paramValNormMean = zeros(size(paramValNorm,1),size(paramValNorm,2));
    paramValNormStd = zeros(size(paramValNorm,1),size(paramValNorm,2));
    paramValNormSEM = zeros(size(paramValNorm,1),size(paramValNorm,2));
    % loop through all time points
    for i = 1:size(paramValNorm,1)
        % loop through all legs
        for j = 1:size(paramValNorm,2)
            thisTandLegs = paramValNorm(i,j,:);

            % remove NaNs
            thisTandLegs(isnan(thisTandLegs)) = [];
            
            paramValNormMean(i,j) = mean(thisTandLegs);
            paramValNormStd(i,j) = std(thisTandLegs);
            paramValNormSEM(i,j) = std(thisTandLegs) / ...
                sqrt(length(thisTandLegs));
        end
    end

    if (semError)
        paramErr = paramValNormSEM;
    else
        paramErr = paramValNormStd;
    end

    % get cell array for violin plot for fwhm
    % preallocate
    fwhmCell = {};
    for i = 1:size(fwhmAll,1)
        fwhmCell{i} = fwhmAll(i,~isnan(fwhmAll(i,:)));
    end

    % plot normalized step param, ipsi and contra legs onto of each other
    % initialize figure
    figure;
    c = colormap('lines');

    for i = 1:3
        subplot(3,1,i);
        hold on;

        % plot
%         errorbar(legT', paramMean(:,i), paramErr(:,i), ...
%             'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(1,:));

        thisNormMean = (paramValNormMean(:,i) - min(paramValNormMean(:,i))) / ...
            (max(paramValNormMean(:,i)) - min(paramValNormMean(:,i)));
        thisNormErr = (paramErr(:,i)) / ...
            (max(paramValNormMean(:,i)) - min(paramValNormMean(:,i)));

        plot_err_patch_v2(legT',thisNormMean,thisNormErr,c(1,:) * 0.8,c(1,:));

%         errorbar(legT', paramMean(:,i+3), paramErr(:,i+3), ...
%             'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(2,:));

        thisNormMean = (paramValNormMean(:,i+3) - min(paramValNormMean(:,i+3))) / ...
            (max(paramValNormMean(:,i+3)) - min(paramValNormMean(:,i+3)));
        thisNormErr = (paramErr(:,i+3)) / ...
            (max(paramValNormMean(:,i+3)) - min(paramValNormMean(:,i+3)));

        plot_err_patch_v2(legT',thisNormMean,thisNormErr,c(2,:) * 0.8,c(2,:));


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
        legend({'','Ipsi', '','Contra'});
    end

    sgtitle(whichParam);



    % plot violin plots of fwhm
    figure;
    for i = 1:3
        subplot(3,1,i);
        hold on;

        violinInd = [i, i+3];

        violin(fwhmCell(violinInd), 'facecolor', c(1:2,:), 'plotlegend',0,...
            'medc',[]);
        xticks(1:2);
        xticklabels({'Ipsi','Contra'});

        % label legs
        if (i == 1)
            title('Front Legs');
        elseif (i == 2)
            title('Mid legs');
        else
            title('Hind legs');
        end

        ylim(yScale);
    end

    sgtitle(whichParam);


    % paired t-test for each set of legs
    for i = 1:3
        thisDiff = fwhmAll(i,:) - fwhmAll(i+3,:);
        [~,p(i)] = ttest(thisDiff);
    end

end