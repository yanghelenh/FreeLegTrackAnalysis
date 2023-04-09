% preprocessBigRigFlyData.m
%
% Script to load big-rig data from Luke and Ashley and output processed
%  data files (one for each fly)
% Also, calls ffmpeg to adjust levels in video files, convert to mp4, 
%  rename to include trial name, and put all videos in one folder
% Processed data contains:
%
% Adapted from main-for-helen jupyter notebook from Luke
% Foraging_for_Helen jupyter notebook from Ashley
% Modified from loadBigRigFlyData.m
%
% Note: ffmpeg command hard coded in only works on mac (codec is named
%  differently on windows)
% 
% CREATED: 3/28/23 - HHY
% UPDATED: 3/28/23 - HHY

%% paths
legVidPath = '/Users/hyang/Dropbox (HMS)/FreeWalkLegAnalysis-Helen/mp4Vids';
pDataPath = '/Users/hyang/Dropbox (HMS)/FreeWalkLegAnalysis-Helen/pData';

% path to ffmpeg folder
binPath = '/usr/local/bin/';

%% define data paths - Luke's data

exptPath = '/Users/hyang/Dropbox (Personal)/Clandinin Lab/Luke_bigrig';

% expts = {'exp-20181102-165424', ...
%          'exp-20181102-175232', ...
%          'exp-20181103-184106', ...
%          'exp-20181104-162518', ...
%          'exp-20181105-115608', ...
%          'exp-20181107-181316', ...
%          'exp-20181108-111101', ...
%          'exp-20181108-143044', ...
%          'exp-20181109-084314'};

expts = {...
         'exp-20181107-181316', ...
         'exp-20181108-111101', ...
         'exp-20181108-143044', ...
         'exp-20181109-084314'};

ffmpegLevels = '"curves=all=''0/0 0.75/1 1/1''"';

curDir = pwd;

%% definte data paths - Ashley's data

exptPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/FreeWalkLegTrack/FromAshley';

expts = {'exp-20220805-101320', 'exp-20220808-114232'};

ffmpegLevels = '"curves=all=''0/0 0.25/1 1/1''"';

curDir = pwd;
 
%% load data

% loop through all data
for i = 1:length(expts)
    % all files and folders in experiment folder
    trialFolders = dir([exptPath filesep expts{i} filesep 'trial*']);
    
    % loop through all trial folders in this experiment folder
    for j = 1:length(trialFolders)
        cd([trialFolders(j).folder filesep trialFolders(j).name]);
        
        % trial name - is name of trial folder
        thisTrialName = trialFolders(j).name;

        % get text files associated with this trial
        trialTxtFiles = dir('*.txt');

        % check that this trial has cnc.txt, cam.txt, and cam_compr.mkv
        if (isfile('cam.txt') && isfile('cnc.txt') && isfile('cam_compr.mkv'))
            % read in cnc data
            cncFile = fopen('cnc.txt');
            cncDat = textscan(cncFile, '%f %f %f\n', ...
                'Delimiter', ',', 'HeaderLines', 1);
            fclose(cncFile);

            trialDat.cnc.t = cncDat{1};
            trialDat.cnc.x = cncDat{2};
            trialDat.cnc.y = cncDat{3};

            % read in cam data
            camFile = fopen('cam.txt');
            camDat = textscan(camFile, '%f %s %f %f %f %f %f\n', ...
                'Delimiter', ',', 'HeaderLines', 1);
            fclose(camFile);

            trialDat.cam.t = camDat{1};
            trialDat.cam.flyPresent = strcmpi(camDat{2}, 'True');
            trialDat.cam.x = camDat{3};
            trialDat.cam.y = camDat{4};
            trialDat.cam.ang = camDat{7};

            % read in stimulus data if it exists
            if(isfile('stimuli.txt'))
                % read in stimulus data
                stimFile = fopen('stimuli.txt');
                stimTemp = textscan(stimFile, '%s', 'Delimiter', {':',','});
                for k = 3:2:(length(stimTemp{1})-2)
                    % extract variable name
                    varName = stimTemp{1}{k}(2:(end-1));
                    % extract variable value
                    if (contains(stimTemp{1}{k+1}, '"'))
                        varVal = stimTemp{1}{k+1}(2:(end-1));
                    else
                        varVal = str2num(stimTemp{1}{k+1});
                    end
    
                    stimDat.(varName) = varVal;
                end
                trialDat.stim = stimDat;
    
                fclose(stimFile);
            % empty vector for trialDat.stim if no stimulus data exists    
            else
                trialDat.stim = [];
            end

            % read in opto data if it exists - get time led on, led off
            if(isfile('opto.txt'))
                % read in file line by line
                optoLines = readlines('opto.txt');

                % initialize opto start and end times
                optoStartT = [];
                optoEndT = [];

                % loop through all lines, find ones with LED info
                for k = 1:length(optoLines)
                    % check if this is a line with LED info
                    if (contains(optoLines(k), 'led'))
                        % split this line at commas
                        thisLineStrs = split(optoLines(k), ', ');

                        % get time (2nd value) - convert to double
                        thisTime = double(thisLineStrs(2));

                        % check whether this is led on or off
                        if (strcmpi(thisLineStrs(3), 'on'))
                            optoStartT = [optoStartT; thisTime];
                        elseif (strcmpi(thisLineStrs(3), 'off'))
                            optoEndT = [optoEndT; thisTime];
                        end
                    end
                end
                trialDat.opto.startT = optoStartT;
                trialDat.opto.endT = optoEndT;
            % empty vector for trialDat.opto if no opto data exists    
            else
                trialDat.opto = [];
            end 

            % call ffmpeg: adjust levels on video, convert to mp4, name
            %  trial name

            % get full output video path
            vidPath = ['''' legVidPath filesep thisTrialName '.mp4'''];

            % generate command for creating ffmpeg video 
            vidCmd = ...
                sprintf('ffmpeg -i cam_compr.mkv -vf %s -pix_fmt yuv420p -b:v 10000k -c:v h264_videotoolbox %s', ...
                ffmpegLevels, vidPath);
            % run command for creating ffmpeg video
            createVidStatus = system(['export PATH=' binPath ' ; ' vidCmd]);
        
            % display whether video file generated successfully
            if ~(createVidStatus)
                fprintf('Video for %s created successfully! \n', ...
                    thisTrialName);
            else
                fprintf('Error creating video for %s. \n', thisTrialName);
            end

            % save trialDat into pData file
            pDataFullPath = [pDataPath filesep thisTrialName '.mat'];
            save(pDataFullPath, 'trialDat','-v7.3');
        end
    end
end
