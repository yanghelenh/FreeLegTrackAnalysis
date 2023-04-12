% extractStepsFromPData.m
%
% Function that computes leg step parameters from raw leg tracking data 
%  (output of loadFreeWalkTrk2PData())
% User selects one or more pData files through GUI
% The pData file must contain the moveNotMove struct (output of
%  extractMoveNotMoveFromPData() and the legTrack struct 
%  (output of loadFreeWalkTrk2PData())
% Saves back into same pData file new structs legSteps, stanceStepParams,
%  and swingStepParams
%
% INPUTS:
%   pDataPath - full path to pData folder
%   peakParams - struct of parameters to use for peak calling, if empty,
%       use default parameters
%
% OUTPUTS: none, but save legSteps, stanceStepParams, and swingStepParams
%   back into pData file
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%
function extractStepsFromPData(pDataPath, peakParams)

end