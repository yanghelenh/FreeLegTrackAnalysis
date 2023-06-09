% plotCondIndptStepParam_allFiles.m
%
% Function to plot 1D step parameter (i.e. not AEP or PEP) value across
%  conditions. Each condition mean +/- SEM or std. Separate subplots
%  for each leg. Specify one step parameter to plot. Select data (output 
%  of saveLegStepParamCond_indpt()) through GUI
%
% INPUTS:
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 6/9/23 - HHY
%
% UPDATED:
%   6/9/23 - HHY
%
function plotCondIndptStepParam_allFiles(whichParam, swingOrStance, ...
    datDir, semError, yScale)
    
    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 
    numLegs = length(subInd);

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
    allMeans = zeros(numFiles, numLegs);
    allErrs = zeros(numFiles, numLegs);
    numBouts = zeros(numFiles, numLegs);
    legendStr = cell(numFiles, numLegs);

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

        % load data
        if(swingOrStance) % stance
            load(outputFullPath, 'selStanceParams', 'legIDs', 'cond');
            paramVal = selStanceParams.(whichParam);
        else % swing
            load(outputFullPath, 'selSwingParams', 'legIDs', 'cond');
            paramVal = selSwingParams.(whichParam);
        end

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
                paramValRad = deg2rad(wrapTo180(rmoutliers(...
                    wrapTo360(paramVal(thisLegLog)))));
                % get mean
                if(swingOrStance) % if stance, wrap to 360 instead, for plotting
                    allMeans(i,j) = wrapTo360(rad2deg(...
                        circ_mean(paramValRad)));
                else % swing, don't wrap
                    allMeans(i,j) = rad2deg(...
                        circ_mean(paramValRad));
                end
                % get std
                thisStd = rad2deg(circ_std(paramValRad));
                % get n
                thisN = length(paramValRad);
            else
                allMeans(i,j) = mean(rmoutliers(paramVal(thisLegLog)));
                thisStd = std(rmoutliers(paramVal(thisLegLog)));
                thisN = length(rmoutliers(paramVal(thisLegLog)));
            end
            
            numBouts(i,j) = thisN;

            % compute errors
            if (semError) % SEM
                allErrs(i,j) = thisStd / sqrt(thisN);
            else % std dev
                allErrs(i,j) = thisStd;
            end

            % get legend str, named based on condition
            thisLegendStr = [];
            for k = 1:length(cond.whichParam)
                if (k>1) && (strcmpi(cond.whichParam{k},cond.whichParam{k-1}))
                    thisLegendStr = [thisLegendStr cond.cond{k}];
                else
                    thisLegendStr = [thisLegendStr cond.whichParam{k}...
                        cond.cond{k}];
                end
                if (k~=length(cond.whichParam))
                    thisLegendStr = [thisLegendStr ', '];
                end
            end

            % add n to str
            thisLegendStr = [thisLegendStr '; n=' num2str(thisN)];

            legendStr{i,j} = thisLegendStr;
        end
    end

    % x indices for plotting
    xInd = 1:numFiles;

    % initialize figure
    f = figure;

    for i = 1:6
        subplot(3,2,subInd(i));
        hold on;

        % plot 
        errorbar(xInd, allMeans(:,i), allErrs(:,i), '.', 'Marker', '_', ...
            'LineWidth', 2);

        % label each condition on x-axis
        set(gca, 'xTick', xInd', 'XTickLabel', legendStr(:,i));

        % axis scale and label
        ylim(yScale);
        xlim([xInd(1)-1, xInd(end)+1]);
        ylabel(whichParam);
    end

    sgtitle(whichParam);
end