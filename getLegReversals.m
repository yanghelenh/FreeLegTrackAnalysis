% getLegReversals.m
%
% Function for getting max and min in leg positions, when leg reverses
%  direction
% Adapted from nonInteractGetLegReversals(), but uses findpeaks to get leg
%  reversals instead of findLegReversals()
% 
% INPUTS:
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   moveNotMove - struct of moving/not-moving bout indicies, starts/ends
%   legRevParams - struct of parameters for getting these leg reversals
%       minProm - parameter MinPeakProminence of findpeaks
%       minDist - parameter MinPeakDistance of findpeaks
%   legIDs - struct of parameters, IDs legs
%       ind - indicies of legs, into raw position matricies, carried
%           throughout
%       names - names of each leg, matching with ind
%
% OUTPUTS:
%   maxIndsAll - indicies (frames) of leg position maxima, all legs in 1
%       vector
%   minIndsAll - indicies of leg position minima
%   maxWhichLeg - number indicating which leg max indicies belong to, same
%       size & matched with maxIndsAll
%   minWhichLeg - which legs for min indicies
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%   6/21/23 - HHY - break out main part of function into helper function 
% 
function [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
    getLegReversals(legTrack, moveNotMove, legRevParams, legIDs)

    % initialize maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg
    maxIndsAll = [];
    minIndsAll = [];
    maxWhichLeg = [];
    minWhichLeg = [];

    % loop through all legs
    for i = 1:length(legIDs.ind)
        thisLegInd = legIDs.ind(i);

        legXPos = legTrack.srnfLegX(:,thisLegInd);

        [maxInds, minInds] = getLegReversals_1Leg(legXPos, ...
            moveNotMove.notMoveInd, legRevParams);

        % return values for this, concatenate into vector for all legs
        maxIndsAll = [maxIndsAll; maxInds];
        minIndsAll = [minIndsAll; minInds];
        
        % which leg marker, generate vectors
        thisMaxWhichLeg = ones(size(maxInds)) * legIDs.ind(i);
        thisMinWhichLeg = ones(size(minInds)) * legIDs.ind(i);
        % concatenate into vector for all legs
        maxWhichLeg = [maxWhichLeg; thisMaxWhichLeg];
        minWhichLeg = [minWhichLeg; thisMinWhichLeg];
    end
end