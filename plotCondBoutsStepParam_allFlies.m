% plotCondBoutsStepParam_allFlies.m
%
% Function to plot 1D step parameters for bouts aligned to yaw velocity
%  peaks. From output of saveLegStepParamCond_bouts()
% One figure with 6 subplots, 1 for each leg
% Can select more than 1 output file, will plot on same graph
% Has option to plot reference steps (from saveLegStepParamCond_indpt()).
%  Select one file through GUI. Plots mean and std/SEM (as specified by
%  input)
% Select output files through GUI
% Note: make sure all output files selected have same maxNumSteps. Plotting
%  will get confused otherwise
%
% INPUTS:
%   whichParam - name of parameter to plot
%   swingOrStance - boolean for whether to plot param during swing (false)
%       or stance (true)
%   datDir - directory with output files
%   yScale - y scale for plots, as [min max]
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
%   6/5/23 - HHY add reference steps
%   6/6/23 - HHY - for reference file, make sure stepDirections mean and
%       std use circular statistics
%   6/7/23 - HHY - remove outliers from reference file
%
function plotCondBoutsStepParam_allFlies(whichParam, swingOrStance, ...
    datDir, yScale, semError, ref)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % all step parameters that are circular variables - need to use
    %  circular stats - 6/6/23 - HHY
    circStepParams = {'stepDirections'};

    % prompt user to select output files from saveLegStepParamCond_bouts()
    disp('Select output files from saveLegStepParamCond_bouts()');
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
        %  saveLegStepParamCond_indpt()
        disp('Select reference file');
        [refOutputName, refOutputPath] = uigetfile('*.mat', ...
            'Select reference file', datDir, 'MultiSelect', 'off');

        % load reference
        refOutputFullPath = [refOutputPath filesep refOutputName];

        if(swingOrStance) % stance
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

        % load data - with appropriate error and phase
        if (semError && swingOrStance) % SEM and Stance
            load(outputFullPath, 'stanceParamMeans', 'stanceParamSEM', ...
                'stanceParamN', 'maxNumSteps', 'cond');

            % if circular parameter and stance, wrap to 360 for plotting
            if(any(strcmpi(whichParam, circStepParams)))
                allMeans{i} = wrapTo360(stanceParamMeans.(whichParam));
            else % not circular
                allMeans{i} = stanceParamMeans.(whichParam);
            end

            allErrs{i} = stanceParamSEM.(whichParam);
            numBouts{i} = stanceParamN.(whichParam);
        elseif (semError && ~swingOrStance) % SEM and Swing
            load(outputFullPath, 'swingParamMeans', 'swingParamSEM', ...
                'swingParamN', 'maxNumSteps', 'cond');

            allMeans{i} = swingParamMeans.(whichParam);
            allErrs{i} = swingParamSEM.(whichParam);
            numBouts{i} = swingParamN.(whichParam);
        elseif (~semError && swingOrStance) % std dev and stance
            load(outputFullPath, 'stanceParamMeans', 'stanceParamStd', ...
                'stanceParamN', 'maxNumSteps', 'cond');

            % if circular parameter and stance, wrap to 360 for plotting
            if(any(strcmpi(whichParam, circStepParams)))
                allMeans{i} = wrapTo360(stanceParamMeans.(whichParam));
            else % not circular
                allMeans{i} = stanceParamMeans.(whichParam);
            end

            allErrs{i} = stanceParamStd.(whichParam);
            numBouts{i} = stanceParamN.(whichParam);
        elseif (~semError && ~swingOrStance) % std dev and swing
            load(outputFullPath, 'swingParamMeans', 'swingParamStd', ...
                'swingParamN', 'maxNumSteps', 'cond');

            allMeans{i} = swingParamMeans.(whichParam);
            allErrs{i} = swingParamStd.(whichParam);
            numBouts{i} = swingParamN.(whichParam);
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
        legendStr{i} = thisLegendStr;
    end

    % if reference, compute means and errors
    if (ref)
        % loop through each leg
        for j = 1:numLegs
            thisLeg = legIDs.ind(j);
            if(swingOrStance)
                thisLegLog = selStanceParams.stepWhichLeg == thisLeg;
            else
                thisLegLog = selSwingParams.stepWhichLeg == thisLeg;
            end
            
            % compute mean, std; check if circular param
            if(any(strcmpi(whichParam, circStepParams)))
                % convert to radians
                refParamValRad = deg2rad(wrapTo180(rmoutliers(wrapTo360(refParamVal(thisLegLog)))));
                % get mean
                if(swingOrStance) % if stance, wrap to 360 instead, for plotting
                    refMeans(j) = wrapTo360(rad2deg(...
                        circ_mean(refParamValRad)));
                else % swing, don't wrap
                    refMeans(j) = rad2deg(...
                        circ_mean(refParamValRad));
                end
                % get std
                thisStd = rad2deg(circ_std(refParamValRad));
                % get n
                thisN = length(refParamValRad);
            else
                refMeans(j) = mean(rmoutliers(refParamVal(thisLegLog)));
                thisStd = std(rmoutliers(refParamVal(thisLegLog)));
                thisN = length(rmoutliers(refParamVal(thisLegLog)));
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

    % get step time points (for x-axis): 0 for step at peak
    stepTPts = [fliplr((1:maxNumSteps) * -1), 0, 1:maxNumSteps];

    % if reference, generate lines for plotting reference values across all
    %  time points
    if(ref)
        refMeansRep = repmat(refMeans,1,length(stepTPts));
        refErrsPosRep = repmat(refMeans + refErrs,1,length(stepTPts));
        refErrsNegRep = repmat(refMeans - refErrs,1,length(stepTPts));
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
            plot(stepTPts, refMeansRep(i,:), 'LineWidth',1, ...
                'Color', 'black');
            % plot error lines
            plot(stepTPts, refErrsPosRep(i,:), 'LineWidth',1, ...
                'Color', 'black', 'LineStyle',':');
            plot(stepTPts, refErrsNegRep(i,:), 'LineWidth',1, ...
                'Color', 'black', 'LineStyle',':');

            % legend
            thisLegendStr = ['Reference; n = ' num2str(refN(i))];
    
            legendStrLegs{1,i} = thisLegendStr;
            % omit legend for error lines
            legendStrLegs{2,i} = '';
            legendStrLegs{3,i} = '';
    
            hold on;
        end

        % loop through all conditions
        for j = 1:numFiles
            % plot for this condition
            errorbar(stepTPts', allMeans{j}(:,i), allErrs{j}(:,i), ...
                'Marker', 'x','LineWidth',1, 'CapSize', 0, 'Color', c(j,:));

            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j} '; n = ' ...
                num2str(numBouts{j}(:,i)')];

            % legend depends on if there's reference
            if(ref)
                legendStrLegs{j+3,i} = thisLegendStr;
            else
                legendStrLegs{j,i} = thisLegendStr;
            end
        end

        % axis scale and label
%         ylim(yScale);
        xScale = xlim;
        xScale(1) = xScale(1) - (0.1 * (stepTPts(end)-stepTPts(1)));
        xScale(2) = xScale(2) + (0.1 * (stepTPts(end)-stepTPts(1)));
        xlim(xScale);

        xticks(stepTPts);

        xlabel('Steps');
        ylabel(whichParam);

        % legend
        legend(legendStrLegs{:,i});
    end

    sgtitle(whichParam);
end