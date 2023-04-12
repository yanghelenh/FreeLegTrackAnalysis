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
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%
function extractMoveNotMoveFromPData(pDataPath)

end