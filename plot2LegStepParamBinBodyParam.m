% plot2LegStepParamBinBodyParam.m
%
% Function to plot two leg step parameters against each other, while
%  binning by a third parameter. Operates on output of
%  saveLegStepParamCond_indpt().
% Assumes saveLegStepParamCond_indpt() output file has values across the
%  range of the bodytraj parameter (i.e. don't use files pre-filtered for 
%  specific values.
% Select saveLegStepParamCond_indpt() output file through GUI
% NOTE: 9/5/23 - don't use step directions as one of the X or Y parameters
%  (currently, this function doesn't do circular stats)
%
% INPUTS:
%   stepParamX - name of leg step parameter that will be on x-axis
%   stepParamY - name of leg step parameter that will be on y-axis
%   whichPhaseX - leg phase for stepParamX
%   whichPhaseY - leg phase for stepParamY
%   binParam - name of parameter to bin by
%   binParamPhase - leg phase for binParam
%   binRange - 2 element vector for [start, end] of range of bodyParam to
%       bin over
%   numBins - number of bins of bodyParam
%   xScale - x-axis scale
%   yScale - y-axis scale
%   datDir - folder with saveLegStepParamCond_indpt() output file
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/5/23 - HHY
%
% UPDATED:
%   9/5/23 - HHY
%
function plot2LegStepParamBinBodyParam(stepParamX, stepParamY, ...
    whichPhaseX, whichPhaseY, binParam, binParamPhase, binRange, ...
    numBins, datDir)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % get bins
    binEdges = linspace(binRange(1), binRange(2), numBins + 1);
    binStarts = binEdges(1:(end-1));
    binEnds = binEdges(2:end);
    binMids = (binStarts + binEnds) / 2;

    % preallocate (6 for number of legs)
    % counts for number of steps in each bin
    binCounts = zeros(numBins, 6);
    % sums for summed value of X/Y parameter for that bin
    binValsX = cell(numBins, 6);
    binValsY = cell(numBins, 6);

     % prompt user to select output file from saveLegStepParamCond_indpt()
    [outputFName, outputPath] = uigetfile('*.mat', ...
        'Select saveLegStepParamCond_indpt() file', ...
        datDir, 'MultiSelect', 'off');

    fullFilePath = [outputPath filesep outputFName];

    % load data
    load(fullFilePath, 'selStanceParams', 'selSwingParams', 'legIDs');

    % get parameters for X, Y, binning
    if (strcmpi(whichPhaseX,'stance'))
        xParamVals = selStanceParams.(stepParamX);
    elseif (strcmpi(whichPhaseX, 'swing'))
        xParamVals = selSwingParams.(stepParamX);
    end

    if (strcmpi(whichPhaseY,'stance'))
        yParamVals = selStanceParams.(stepParamY);
    elseif (strcmpi(whichPhaseY, 'swing'))
        yParamVals = selSwingParams.(stepParamY);
    end

    if (strcmpi(binParamPhase,'stance'))
        binParamVals = selStanceParams.(binParam);
    elseif (strcmpi(binParamPhase, 'swing'))
        binParamVals = selSwingParams.(binParam);
    end


    % loop through all steps, put into appropriate bin
    for i = 1:length(binParamVals)
        thisLegInd = selStanceParams.stepWhichLeg(i);

        % get which bin this belongs to
        thisBinInd = find((binParamVals(i) >= binStarts) & ...
            (binParamVals(i) < binEnds));

        % update bin counts and vals appropriately
        if ~isempty(thisBinInd)
            binCounts(thisBinInd,thisLegInd) = ...
                binCounts(thisBinInd,thisLegInd) + 1;

            binValsX{thisBinInd,thisLegInd} = ...
                [binValsX{thisBinInd,thisLegInd}; xParamVals(i)];

            binValsY{thisBinInd,thisLegInd} = ...
                [binValsY{thisBinInd,thisLegInd}; yParamVals(i)];
        end
    end

    % get mean and SEM for each bin
    binMeanX = zeros(numBins,6);
    binMeanY = zeros(numBins,6);
    binSEMx = zeros(numBins,6);
    binSEMy = zeros(numBins,6);

    % loop through all bins
    for i = 1:numBins
        % loop through all legs
        for j = 1:6
            binMeanX(i,j) = mean(rmoutliers(binValsX{i,j}));
            binMeanY(i,j) = mean(rmoutliers(binValsY{i,j}));
            binSEMx(i,j) = std(rmoutliers(binValsX{i,j})) / ...
                sqrt(length(rmoutliers(binValsX{i,j})));
            binSEMy(i,j) = std(rmoutliers(binValsY{i,j})) / ...
                sqrt(length(rmoutliers(binValsY{i,j})));
        end
    end

    % plot
    figure;
    c = colormap('cool');

    % get colormap indices, for numBins
    colorSpace = floor(size(c,1) / numBins);

    colorInd = 1:colorSpace:size(c,1);

    
    for i = 1:6
        
        subplot(3,2,subInd(i));
        hold on;

        % plot one bin at a time so colors can be different
        for j = 1:numBins
            % plot errorbar
            errorbar(binMeanX(j,i), binMeanY(j,i), binSEMy(j,i), ...
                binSEMy(j,i), binSEMx(j,i), binSEMx(j,i), ...
                'Color', c(colorInd(j),:));
        end

%         xlim(xScale);
%         ylim(yScale);

        % label x axis
        if (i == 3 || i == 6)
            xlabel(stepParamX);
        end

        % label y axis
        if (i == 4 || i == 5 || i == 6)
            ylabel(stepParamY);
        end

        title(legIDs.name{i});
    end

    sgtitle(binParam)
end

