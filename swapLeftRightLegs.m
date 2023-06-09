% swapLeftRightLegs.m
%
% Helper function to swap left and right leg indices
%
% INPUTS:
%   legInd - current index of leg
%   matchedLegInd - matrix matching legs, one column for left, one column
%       for right; same row is same leg
%
% OUTPUTS:
%   newLegInd - new leg index, swapping left and right
%
% CREATED: 6/9/23 - HHY
%   
% UPDATED: 
%   6/9/23 - HHY
%
function newLegInd = swapLeftRightLegs(legInd, matchedLegInd)
    % flip left and right legs
    % r for row index, c for column index
    [r, c] = ind2sub(size(matchedLegInd),...
        find(matchedLegInd == legInd));
    % matched ind are across columns, flip columns
    if (c==1)
        c = 2;
    else
        c = 1;
    end
    % new leg index, swapping left and right
    newLegInd = matchedLegInd(r,c);

end