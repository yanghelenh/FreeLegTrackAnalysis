% plotContLegStepParamCond_bouts.m
%
% Function to plot output of saveContLegStepParamCond_bouts(), compared
%  with reference from saveContLegStepParamCond_indpt()
% One figure with 6 subplots, 1 for each leg
% Can select more than 1 output file, will plot on same graph
% Has option to plot reference steps (from saveLegStepParamCond_indpt()).
%  Select one file through GUI. Plots mean and std/SEM (as specified by
%  input)
% Select output files through GUI
% Note: make sure all output files selected have same legT. Plotting
%  will get confused otherwise
%
% INPUTS:
%   whichParam - name of parameter to plot
%   datDir - directory with output files
%   yScale - y scale for plots, as [min max]
%   semError - boolean for whether to plot SEM; if false, plots std dev
%   ref - boolean for whether to plot reference
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/6/23 - HHY
%
% UPDATED:
%   9/6/23 - HHY
%
function plotContLegStepParamCond_bouts(whichParam, datDir, yScale, ...
    semError, ref)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % all step parameters that are circular variables - need to use
    %  circular stats
    circStepParams = {'stepDirections'};

    % prompt user to select output files from saveContLegStepParamCond_bouts()
    disp('Select output files from saveContLegStepParamCond_bouts()');
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select cond_bouts files', datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFiles = length(outputFNames);
    else
        numFiles = 1;
    end

    % preallocate
    allMeans = cell(numFiles, 1);
    allErrs = cell(numFiles, 1);
    numBouts = cell(numFiles, 1);
    legendStr = cell(numFiles, 1);

    % legend depends on if there's reference (+3 for ref line and error
    %  lines)
    if(ref)
        legendStrLegs = cell(numFiles+3, 6);
    else
        legendStrLegs = cell(numFiles, 6);
    end

    % if reference
    if (ref)
        % prompt user to select reference file from 
        %  saveContLegStepParamCond_indpt()
        disp('Select reference file');
        [refOutputName, refOutputPath] = uigetfile('*.mat', ...
            'Select reference file', datDir, 'MultiSelect', 'off');

        % load reference
        refOutputFullPath = [refOutputPath filesep refOutputName];

        load(refOutputFullPath, 'selLegStepsCont', 'legIDs');
        
        refParamVal = selLegStepsCont.(whichParam);

        % number of legs
        numLegs = length(legIDs.ind);

        % preallocate
        refMeans = zeros(numLegs, 1);
        
        refErrs = zeros(numLegs, 1);

        refN = zeros(1, numLegs);
    end

    % loop through all files, save output matrices for each file
    % save also, cond, numBouts
    for i = 1:numFiles
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data - with appropriate error
        if (semError) % SEM
            load(outputFullPath, 'selLegStepsCont', 'legStepsContMean', ...
                'legStepsContSEM', 'legT', 'cond', 'fwdVelCond');

            % if circular parameter, wrap to 360 for plotting
            if(any(strcmpi(whichParam, circStepParams)))
                allMeans{i} = wrapTo360(legStepsContMean.(whichParam));
            else % not circular
                allMeans{i} = legStepsContMean.(whichParam);
            end

            allErrs{i} = legStepsContSEM.(whichParam);
            numBouts{i} = size(selLegStepsCont.(whichParam), 1);
        elseif (~semError) % std dev
            load(outputFullPath, 'selLegStepsCont', 'legStepsContMean', ...
                'legStepsContStd', 'legT', 'cond', 'fwdVelCond');

            % if circular parameter, wrap to 360 for plotting
            if(any(strcmpi(whichParam, circStepParams)))
                allMeans{i} = wrapTo360(legStepsContMean.(whichParam));
            else % not circular
                allMeans{i} = legStepsContMean.(whichParam);
            end

            allErrs{i} = legStepsContStd.(whichParam);
            numBouts{i} = size(selLegStepsCont.(whichParam), 3);
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
        % add fwdVelCond, change in forward vel, to legend string
        fwdCondStr = sprintf('\\Delta Fwd Vel %.3f-%.3f',...
            fwdVelCond.change(1), fwdVelCond.change(2));
        thisLegendStr = [thisLegendStr ', ' fwdCondStr];
        legendStr{i} = thisLegendStr;
    end

    % if reference, compute means and errors
    if (ref)
        % loop through each leg
        for j = 1:numLegs
            thisRefParamVal = refParamVal(:,j);
            thisRefParamVal(isnan(thisRefParamVal)) = [];
            
            % compute mean, std; check if circular param
            if(any(strcmpi(whichParam, circStepParams)))
                % convert to radians
                refParamValRad = deg2rad(thisRefParamVal);
                % get mean
                refMeans(j) = wrapTo360(rad2deg(...
                    circ_mean(refParamValRad)));
                % get std
                thisStd = rad2deg(circ_std(refParamValRad));
                % get n
                thisN = length(refParamValRad);
            else
                refMeans(j) = mean(thisRefParamVal);
                thisStd = std(thisRefParamVal);
                thisN = length(thisRefParamVal);
            end
            
            refN(j) = thisN;

            % compute errors
            if (semError) % SEM
                refErrs(j) = thisStd / sqrt(thisN);
            else % std dev
                refErrs(j) = thisStd;
            end
    
        end
    end


    % if reference, generate lines for plotting reference values across all
    %  time points
    if(ref)
        refMeansRep = repmat(refMeans,1,length(legT));
        refErrsPosRep = repmat(refMeans + refErrs,1,length(legT));
        refErrsNegRep = repmat(refMeans - refErrs,1,length(legT));
    end

    % initialize figure
    f = figure;
    c = colormap('lines');

    for i = 1:6
        subplot(3,2,subInd(i));
        hold on;

        % if reference, plot as solid line for mean and dotted lines for
        %  error around mean
        if(ref)
            % plot mean line
            plot(legT, refMeansRep(i,:), 'LineWidth',1, ...
                'Color', 'black');
            % plot error lines
            plot(legT, refErrsPosRep(i,:), 'LineWidth',1, ...
                'Color', 'black', 'LineStyle',':');
            plot(legT, refErrsNegRep(i,:), 'LineWidth',1, ...
                'Color', 'black', 'LineStyle',':');

            % legend
            thisLegendStr = ['Reference; n = ' num2str(refN(i))];
    
            legendStrLegs{1,i} = thisLegendStr;

            % omit legend for error lines
            legendStrLegs{2,i} = '';
            legendStrLegs{3,i} = '';
        end

        % loop through all conditions
        for j = 1:numFiles
            % plot for this condition
            errorbar(legT', allMeans{j}(:,i), allErrs{j}(:,i), ...
                'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(j,:));

            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j} '; n = ' ...
                num2str(numBouts{j})];

            % legend depends on if there's reference
            if(ref)
                legendStrLegs{j+3,i} = thisLegendStr;
            else
                legendStrLegs{j,i} = thisLegendStr;
            end
        end

        % axis scale and label

        % for AEP
%         if (subInd(i)==1 || subInd(i)==2)
%             ylim([-0.9 -0.5]);
%         elseif (subInd(i)==3 || subInd(i)==4)
%             ylim([-0.3 0.1]);
%         else
%             ylim([0.15 0.55]);
%         end

        % for PEP
%         if (subInd(i)==1 || subInd(i)==2)
%             ylim([-0.6 -0.2]);
%         elseif (subInd(i)==3 || subInd(i)==4)
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

        % legend
        legend(legendStrLegs{:,i});
    end

    sgtitle(whichParam);
end