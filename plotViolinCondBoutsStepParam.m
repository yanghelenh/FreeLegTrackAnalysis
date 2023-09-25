% plotViolinCondBoutsStepParam.m
%
% Function to plot 1D step parameters for bouts aligned to yaw velocity
%  peaks. From output of saveLegStepParamCond_bouts()
% One figure with 6 subplots, 1 for each leg
% Can select more than 1 output file, will plot on same graph
% Has option to plot reference step (from saveLegStepParamCond_indpt()).
%  Select one file through GUI. Plots mean and std/SEM (as specified by
%  input)
% Note: only works on files with single step. For multiple steps, use
%  plotCondBoutsStepParam_allFlies()
% Select output files through GUI
%
% INPUTS:
%   whichParam - name of parameter to plot
%   swingOrStance - boolean for whether to plot param during swing (false)
%       or stance (true)
%   datDir - directory with output files
%   yScale - y scale for plots, as [min max]
%   ref - boolean for whether to plot reference
%   condNames - cell array of condition names
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/21/23 - HHY
%
% UPDATED:
%   9/21/23 - HHY
%
function plotViolinCondBoutsStepParam(whichParam, swingOrStance, ...
    datDir, yScale, ref, condNames)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % all step parameters that are circular variables - need to use
    %  circular stats - 6/6/23 - HHY
    circStepParams = {'stepDirections'};

    invStepParams = {'stepAEPX', 'stepPEPX'};

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
    allVals = cell(numFiles, 6);

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
            refWhichLeg = selStanceParams.stepWhichLeg;
        else % swing
            load(refOutputFullPath, 'selSwingParams', 'legIDs');
            refParamVal = selSwingParams.(whichParam);
            refWhichLeg = selSwingParams.stepWhichLeg;
        end

        % number of legs
        numLegs = length(legIDs.ind);
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
        if (swingOrStance) % Stance
            load(outputFullPath, 'selStanceParams', ...
                'stanceParamN', 'cond', 'fwdVelCond');

            % get this value
            thisVal = squeeze(selStanceParams.(whichParam));

            numBouts{i} = stanceParamN.(whichParam);
        else
            load(outputFullPath, 'selSwingParams',...
                'swingParamN', 'cond', 'fwdVelCond');

            thisVal = squeeze(swingParamMeans.(whichParam));
            numBouts{i} = swingParamN.(whichParam);

        end

        % loop through all legs, remove outliers and NaNs
        for j = 1:6
            thisValLeg = thisVal(j,:);
            % if circular, and stance, wrap to 360
            if any(strcmpi(whichParam,circStepParams)) && swingOrStance
                thisValLeg = wrapTo360(thisValLeg);
            end

            % remove NaN
            thisValLeg(isnan(thisValLeg)) = [];

            % remove outliers
            thisValLeg = rmoutliers(thisValLeg);

            % add to cell array
            allVals{i,j} = thisValLeg;
        end
    end

    % if reference
    if (ref)
        allRefVals = cell(numLegs,1);
        % loop through each leg
        for j = 1:numLegs
            thisLeg = legIDs.ind(j);
            thisLegLog = refWhichLeg == thisLeg;
            
            % get values
            if(any(strcmpi(whichParam, circStepParams)))
                allRefVals{j} = rmoutliers(wrapTo360(refParamVal(thisLegLog)));
                % get n
                thisN = length(allRefVals{j});
            else
                allRefVals{j} = rmoutliers(refParamVal(thisLegLog));
                thisN = length(allRefVals{j});
            end
            
            refN(j) = thisN;    
        end
    end

    % if reference, compute means
    if (ref)
        % loop through each leg
        for j = 1:numLegs
            thisLeg = legIDs.ind(j);
            thisLegLog = refWhichLeg == thisLeg;
            
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
            else
                refMeans(j) = mean(rmoutliers(refParamVal(thisLegLog)));
            end
            
            refN(j) = thisN; 
        end
    end

    % initialize figure
    f = figure;
    f.Position = [100 100 (length(condNames)+1)*200 800];
    c = colormap('lines');

    for i = 1:6
        violinCell = {};
        subplot(3,2,subInd(i));
        hold on;

        % generate cell arrays for violin plot
        for j = 1:numFiles
            if(ref)
                    if(any(strcmpi(whichParam, circStepParams)))
                        refSubVal = wrapTo180(allVals{j,i} - refMeans(i));
                    else
                        refSubVal = allVals{j,i} - refMeans(i);
                    end
                    violinCell = [violinCell, refSubVal];
    %             for j = 1:numFiles
    %                 violinCell = [violinCell, allVals{j,i}];
    %             end
            else
                violinCell = [violinCell, allVals{j,i}];
            end
        end

        numViolins = length(violinCell);

        % plot violin plots
        violin(violinCell, 'facecolor', c(1:numViolins,:), 'plotlegend',0,...
            'medc',[]);
        xticks(1:numViolins);
        xticklabels(condNames);
        

        % axis scale and label

        % for AEP
%         if (subInd(i)==1 || subInd(i)==2)
%             ylim([-1.1 -0.45]);
%             yticks([-1.1 -0.9 -0.7 -0.5]);
%         elseif (subInd(i)==3 || subInd(i)==4)
%             ylim([-0.5 0.25]);
%             yticks([-0.5 -0.3 -0.1 0.1]);
%         else
%             ylim([-0.1 0.65]);
%             yticks([-0.1 0.1 0.3 0.5]);
%         end

        % for PEP
%         if (subInd(i)==1 || subInd(i)==2)
%             ylim([-0.8 0]);
%         elseif (subInd(i)==3 || subInd(i)==4)
%             ylim([-0.3 0.5]);
%         else
%             ylim([0.15 0.95]);
%         end

        line(xlim,[0,0],'Color','k');

        ylim(yScale);
        ylabel(legIDs.name(i));

        if any(strcmpi(whichParam,invStepParams))
            set(gca,'YDir','reverse');
        end
    end

    sgtitle(whichParam);
end