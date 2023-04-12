% extractMoveNotMoveFromPData.m
%
% Function that extracts moving and not moving time periods (using 
%  body velocities) and then computes leg step parameters from raw leg 
%  tracking data (output of loadFreeWalkTrk2PData())
% User selects single pData file through GUI. This pData file must contain
%  the bodytraj struct (output of preprocessBodyTraj())
% Calls one interactive figure to determine moving/not-moving bouts
%  and the other to correct the leg position max/min calls (needed for
%  determining steps)
%
% INPUTS:
%   pDataPath - full path to pData folder
%   
% OUTPUTS:
%   none, but saves moveNotMove struct into same pData file
%     notMoveInd - indices of all not moving frames
%     notMoveBout - start and end indices of all not moving bouts
%     moveInd - indices of all moving frames
%     moveBout - start and end indices of all moving bouts
%     notMoveParams - final parameter values used
%       transThresh - translational speed threshold, in m/s
%       angThresh - angular speed threshold, in deg/s
%       minBoutLen - minimum bout length, in sec
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%
function extractMoveNotMoveFromPData(pDataPath)

    % prompt user for pData file; defaults to folder containing pData files
    disp('Select pData file to process');
    [pDataName, thisPDataPath] = uigetfile('*.mat', 'Select pData file', ...
        pDataPath);
    % full path to pData file
    pDataFullPath = [thisPDataPath filesep pDataName];

    fprintf('Selected %s\n', pDataName);

    % get variables saved in pData file
    pDatVars = whos('-file', pDataFullPath);

    pDatVarsNames = cell(size(pDatVars));
    
    % convert pDatVars into cell array of just names
    for i = 1:length(pDatVars)
        pDatVarsNames{i} = pDatVars(i).name;
    end

    % check that this pData file has all the variables needed for later
    %  analysis
    if (~any(strcmpi(pDatVarsNames, 'bodytraj')))
        fprintf('%s does not contain the bodytraj struct. Ending \n', ...
            pDataName);
        return;
    end

    % check if this pData file has moveNotMove struct
    % if yes, ask if user wants to rerun
    % end function if user is not rerunning
    if (any(strcmpi(pDatVarsNames, 'moveNotMove')))
        fprintf('Moving/not moving bouts already determined for %s \n', ...
            pDataName);
        % ask if user wants to redo
        contStr = input('Rerun moving/not moving selection? y/n ', 's');
        % not yes
        if ~(strcmpi(contStr, 'y'))
            disp('Not re-running moving/not moving selection. Ending.');
            return;
        end
    end

    % load data
    load(pDataFullPath, 'bodytraj');

    % some parameters
    tRange = 60; % sec, amount of data to display
    xMax = bodytraj.tZeroed(end); % max value

    % initialize defaults for not move parameters
    notMoveParams.transThresh = 0.002; % in m/s
    notMoveParams.angThresh = 20; % in deg/s
    notMoveParams.minBoutLen = 0.5; % in s
    
    % slider parameters
    sldXPos = 1200;
    sldYPosStart = 850;
    sldYPosEnd = 200;
    sldHeight = 20;
    sldWidth = 300;
    numSld = length(fieldnames(notMoveParams));
    sldYSpace = round((sldYPosStart - sldYPosEnd) / numSld);

    % ranges for all notMoveParams
    nmpRanges.transThresh = [0 0.005];
    nmpRanges.angThresh = [0 100];
    nmpRanges.minBoutLen = [0 5];
    
    % initialize logical for done button press
    userDone = 0;
    
    % remember inital parameters
    initNotMoveParams = notMoveParams;


    % interframe interval for bodytraj
    ifi = median(diff(bodytraj.tZeroed));
    sampRate = 1/ifi; % sample rate for bodytraj

    % convert min bout length in sec to samples
    minBoutLenSamp = notMoveParams.minBoutLen * sampRate;

    % translational and angular speed
    transSpdSmo = bodytraj.transVelSmoS;
    angSpdSmo = abs(bodytraj.angVelSmoS);

    % get not moving calls with initial  parameters
    [notMoveInd, notMoveStartInd, notMoveEndInd] = ...
        findFlyNotMovingBody(transSpdSmo, angSpdSmo, ...
        notMoveParams.transThresh, notMoveParams.angThresh, ...
        minBoutLenSamp);

    % get shading for not moving
    notMovingX = [notMoveStartInd'; notMoveStartInd'; ...
        notMoveEndInd'; notMoveEndInd'];
    notMovingXT = bodytraj.tZeroed(notMovingX);

    tY0 = zeros(size(notMoveStartInd'));
    aY0 = zeros(size(notMoveStartInd'));
    tY1 = ones(size(notMoveStartInd')) * 1.1 * max(transSpdSmo);
    aY1 = ones(size(notMoveStartInd')) * 1.1 * max(angSpdSmo);
    tNotMovingY = [tY0; tY1; tY1; tY0];
    aNotMovingY = [aY0; aY1; aY1; aY0];


    % initialize figure
    f = figure('Position', [20 20 1600 920]);

    % plot translational speed
    tAx = subplot('Position', [0.05 0.6 0.6 0.35]);
    plot(tAx,bodytraj.tZeroed, transSpdSmo);
    hold on;
    
    % plot shading for not moving bouts
    patch(tAx,notMovingXT, tNotMovingY, 'black', 'FaceAlpha', 0.3');
    
    xlim([0 tRange]);
    title('Translational Speed (m/s)');
    xlabel('Time (s)');
    
    
    % plot angular speed
    aAx = subplot('Position', [0.05 0.15 0.6 0.35]);
    plot(aAx,bodytraj.tZeroed, angSpdSmo);
    hold on;
    
    % plot shading for not moving bouts
    patch(aAx,notMovingXT, aNotMovingY, 'black', 'FaceAlpha', 0.3');
    
    xlim([0 tRange]);
    title('Angular Speed (deg/s)');
    xlabel('Time (s)');
    
    % link axes for translational and angular speed
%     linkaxes([tAx aAx], 'x');


    % get all notMoveParam names
    nmpNames = fieldnames(notMoveParams);
    % get all slider names
    sldNames = nmpNames;

    % text for names of parameter sliders
    % invisible axes for text
    txtAx = axes('Position',[0 0 1 1], 'Visible', 'off');
    set(gcf, 'CurrentAxes', txtAx);
    for i = 1:length(sldNames)
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        text(sldXPos - 120, thisSldYPos + 10, sldNames{i}, ...
             'Units', 'pixels', 'FontSize', 12);
    end
    % text for values of parameter sliders
    allTxtH = {}; % handles to text obj
    for i = 1:length(sldNames)
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        thisDispVal = num2str(notMoveParams.(sldNames{i}));
        
        allTxtH{i} = text(sldXPos + 320, thisSldYPos + 10, thisDispVal, ...
             'Units', 'pixels', 'FontSize', 12);
    end
    
    % slider adjusting view of leg positions
    tSlider = uicontrol(f, 'Style', 'slider', 'Position', [100 50 600 20]);
    tSlider.Value = 0;
    tSlider.Callback = @updateTLim;
    
    % function to update display region for plot, every time that slider is
    %  moved
    function updateTLim(src, event)
        xlim(tAx, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
        xlim(aAx, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
    end

    % sliders for adjusting each of the parameter values
    % initialize cell array for slider objects
    allSld = {};
    % loop through all parameters, creating sliders
    for i = 1:length(sldNames)
        allSld{i} = uicontrol(f, 'Style', 'slider');
        
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        allSld{i}.Position = [sldXPos thisSldYPos sldWidth sldHeight];
        
        allSld{i}.Value = notMoveParams.(sldNames{i});
        thisParamRange = nmpRanges.(sldNames{i});
        allSld{i}.Max = thisParamRange(2);
        allSld{i}.Min = thisParamRange(1);
        allSld{i}.Callback = {@updateGraph, i};  
    end
    
    % function for updating the figure, notMoveParams, notMove indices
    %  every time a slider is moved
    function updateGraph(src, event, nameInd)

        thisVal = allSld{nameInd}.Value;

        % change the apppropriate parameter value in notMoveParams
        notMoveParams.(sldNames{nameInd}) = thisVal;

        % convert minBoutLen to samples
        minBoutLenSamp = notMoveParams.minBoutLen * sampRate;

        % get not moving calls
        [notMoveInd, notMoveStartInd, notMoveEndInd] = ...
            findFlyNotMovingBody(transSpdSmo, angSpdSmo, ...
            notMoveParams.transThresh, notMoveParams.angThresh, ...
            minBoutLenSamp);
        
        
        % update patches
        % get shading for not moving
        notMovingX = [notMoveStartInd'; notMoveStartInd'; ...
            notMoveEndInd'; notMoveEndInd'];
        notMovingXT = bodytraj.tZeroed(notMovingX);
    
        tY0 = zeros(size(notMoveStartInd'));
        aY0 = zeros(size(notMoveStartInd'));
        tY1 = ones(size(notMoveStartInd')) * 1.1 * max(transSpdSmo);
        aY1 = ones(size(notMoveStartInd')) * 1.1 * max(angSpdSmo);
        tNotMovingY = [tY0; tY1; tY1; tY0];
        aNotMovingY = [aY0; aY1; aY1; aY0];
        
        % plot trans Spd
        cla(tAx);
        plot(tAx,bodytraj.tZeroed, transSpdSmo);
    
        % plot shading for not moving bouts
        patch(tAx,notMovingXT, tNotMovingY, 'black', 'FaceAlpha', 0.3');


        % plot ang Spd
        cla(aAx);
        plot(aAx,bodytraj.tZeroed, angSpdSmo);
    
        % plot shading for not moving bouts
        patch(aAx,notMovingXT, aNotMovingY, 'black', 'FaceAlpha', 0.3');
        
        % update display value around slider
        allTxtH{nameInd}.String = num2str(thisVal);
    end

    % button for user to reset parameter values to original ones
    resetButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Reset', ...
        'Position', [1100 50 50 20]);
    resetButton.Callback = @resetPushed;
    
    % function for resetting notMoveParams when reset button pushed
    function resetPushed(src, event)
        notMoveParams = initNotMoveParams;

        % convert minBoutLen to samples
        minBoutLenSamp = notMoveParams.minBoutLen * sampRate;
        
        % get not moving calls with initial parameters
        [notMoveInd, notMoveStartInd, notMoveEndInd] = ...
            findFlyNotMovingBody(transSpdSmo, angSpdSmo, ...
            notMoveParams.transThresh, notMoveParams.angThresh, ...
            minBoutLenSamp);
        
        % update all slider values
        for j = 1:length(sldNames)
            allSld{j}.Value = notMoveParams.(sldNames{j});
            allTxtH{j}.String = num2str(notMoveParams.(sldNames{j}));
        end

        
        % update patches
        % get shading for not moving
        notMovingX = [notMoveStartInd'; notMoveStartInd'; ...
            notMoveEndInd'; notMoveEndInd'];
        notMovingXT = bodytraj.tZeroed(notMovingX);
    
        tY0 = zeros(size(notMoveStartInd'));
        aY0 = zeros(size(notMoveStartInd'));
        tY1 = ones(size(notMoveStartInd')) * 1.1 * max(transSpdSmo);
        aY1 = ones(size(notMoveStartInd')) * 1.1 * max(angSpdSmo);
        tNotMovingY = [tY0; tY1; tY1; tY0];
        aNotMovingY = [aY0; aY1; aY1; aY0];
        
        % plot trans Spd
        cla(tAx);
        plot(bodytraj.tZeroed, transSpdSmo);
    
        % plot shading for not moving bouts
        patch(notMovingXT, tNotMovingY, 'black', 'FaceAlpha', 0.3');


        % plot ang Spd
        cla(aAx);
        plot(bodytraj.tZeroed, angSpdSmo);
    
        % plot shading for not moving bouts
        patch(notMovingXT, aNotMovingY, 'black', 'FaceAlpha', 0.3');
    end


    % button for when user is done
    doneButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Done', ...
        'Position', [100 30 50 20]);
    doneButton.Callback = @donePushed;
    
    % function for when user is done
    function donePushed(src, event)
        % toggle logical to done, stops loop
        userDone = 1; 
    end
    
    % loop until user hits done button
    while ~userDone
        pause(0.1);
    end


    % get not move as bouts too
    notMoveBout = [notMoveStartInd, notMoveEndInd];

    % get moving ind and bouts
    allBodytrajInd = (1:length(bodytraj.tZeroed))';

    moveInd = allBodytrajInd(~ismember(allBodytrajInd, notMoveInd));
    if ~isempty(moveInd)
        [moveBoutStartInd, moveBoutEndInd, ~] = findBouts(moveInd);
        moveBout = [moveBoutStartInd, moveBoutEndInd];
    else
        moveBout = [];
    end
    
    close(f);

    % save into struct
    moveNotMove.notMoveInd = notMoveInd;
    moveNotMove.notMoveBout = notMoveBout;
    moveNotMove.moveInd = moveInd;
    moveNotMove.moveBout = moveBout;
    moveNotMove.notMoveParams = notMoveParams;

    % update pData file
    save(pDataFullPath, 'moveNotMove', '-append');
end