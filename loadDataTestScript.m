addpath(genpath('/Users/hyang/Documents/EM-analysis-code/'));
bodyCSVPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/FreeWalkLegTrack/melanogaster_legAndBodyCSVs/melanogaster_trial_5_20181011_103827.csv';
legCSVPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/FreeWalkLegTrack/melanogaster_legAndBodyCSVs/position_matrix_melanogaster_trial_5_20181011_103827_cam_compr.csv';
bodyMat = readmatrix(bodyCSVPath, 'NumHeaderLines',1,'Delimiter',',','OutputType','char');
legMat = readmatrix(legCSVPath, 'NumHeaderLines',1','Delimiter',',','OutputType','double');
emptyLeg = isempty(legMat);
emptyLeg = ismissing(legMat);
legMat(5,:)
legMat(6,:)
clearvars
bodyCSVPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/FreeWalkLegTrack/melanogaster_legAndBodyCSVs/melanogaster_trial_5_20181011_103827.csv';
legCSVPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/FreeWalkLegTrack/melanogaster_legAndBodyCSVs/position_matrix_melanogaster_trial_5_20181011_103827_cam_compr.csv';
bodyMat = readmatrix(bodyCSVPath, 'NumHeaderLines',1,'Delimiter',',','OutputType','char');
legMat = readmatrix(legCSVPath, 'NumHeaderLines',1','Delimiter',',','OutputType','double');


%%

combCSVPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/FreeWalkLegTrack/dros_gait_data_for_helen_032223.csv';

combMat = readtable(combCSVPath, 'Delimiter',',', 'ReadVariableNames',true);

melInd = find(strcmpi(combMat.strain,'melanogaster'));

melMat = combMat(melInd,:);

trialNames = unique(melMat.trial_name);
numTrials = length(trialNames);

% separate out contiguous bouts
tDiff = diff(melMat.time);

boutTransInd = find(diff(melMat.time) > 0.02) + 1;

boutTransInd = [1; boutTransInd];

bouts = cell(length(boutTransInd),1);

for i = 1:length(boutTransInd)
    boutStart = boutTransInd(i);
    if (i~=length(boutTransInd))
        boutEnd = boutTransInd(i+1) - 1;
    else
        boutEnd = size(melMat,1);
    end

    bouts{i} = melMat(boutStart:boutEnd,:);
end

testBout = bouts{22};