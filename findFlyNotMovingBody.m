% findFlyNotMovingBody.m
%
% Function that takes in fly's smoothed translational speed and smoothed
%  angular speend, thresholds on those speeds, and a minimum bout length 
%  to extract when the fly is moving or not moving.
% Fly is not moving when speed below both thresholds
% Merges bouts that are too short and assigns 
%
% Helper function for extractMoveNotMoveFromPData()
%
% Adapted from findFlyNotMovingFt()
%
% INPUTS:
%   transSpdSmo - smoothed translational speed
%   angSpdSmo - smoothed angular speed
%   transThresh - threshold on translational speed, below which not moving
%   angThresh - threshold on angular speed, below which not moving
%   minBoutLenSamp - minimum bout length of moving or not moving bout, in
%       samples
%
% OUTPUTS:
%   notMoveInd - logical for each index of transSpdSmo/angSpdSmo, 1 when 
%       fly not moving, 0 when fly moving
%   notMoveStartInd - start indices for when fly not moving
%   notMoveEndInd - end indices for when fly not moving
%
% CREATED: 4/12/23 - HHY
%
% UPDATED:
%   4/12/23 - HHY
%

function [notMoveInd, notMoveStartInd, notMoveEndInd] = ...
    findFlyNotMovingBody(transSpdSmo, angSpdSmo, transThresh, ...
    angThresh, minBoutLenSamp)

    % make row vectors
    if iscolumn(transSpdSmo)
        transSpdSmo = transSpdSmo';
    end

    if iscolumn(angSpdSmo)
        angSpdSmo = angSpdSmo';
    end

    % use threshold to get not moving indices
    % not moving is when translational speed and angular speed are below
    %  their respective thresholds
    notMoveLogical = (transSpdSmo <= transThresh) & ...
        (angSpdSmo <= angThresh);
    
    % indicies where fly transitions between moving and not moving
    transInd = find(diff(notMoveLogical)) + 1;
    
    % add edges to transitions, so first and last bouts are included
    if transInd(end) == length(notMoveLogical)
        transInd = [1 transInd];
    else
        transInd = [1 transInd (length(notMoveLogical)+1)];
    end
    
    % duration of each bout of movement and not movement
    boutDur = diff(transInd);
    
    % deal with bouts that are too short (less than minBoutLen) - merge
    % with other, adjacent short bouts or merge into longer sequence
    whichBout = 1;
    while whichBout<=length(boutDur)
        % bout is too short
        if (boutDur(whichBout) < minBoutLenSamp)
            % index into transInd for start of bout
            boutStartInd = whichBout;
            
            % continue until bouts stop being too short
            for k = (whichBout+1):length(boutDur)
                if (boutDur(k) < minBoutLenSamp)
                    whichBout = k;
                else
                    break;
                end
            end
            
            % index into notMoveLogical for bout transitions
            boutStartMLInd = transInd(boutStartInd);
            % index into transInd for end of bout
            boutEndInd = whichBout +1;
            % equivalent for index into notMoveLogical
            boutEndMLInd = transInd(boutEndInd) - 1;
            
            % is this a moving or not moving bout
            boutMoveLogical = notMoveLogical(boutStartMLInd);
            
            % multiple short bouts, to be merged
            if (whichBout ~= boutStartInd)
                % assign all of these short bouts to 1 longer bout, of the
                %  predominant type (more samples)
                % get number of moving and not moving samples during this
                %  window
                numNotMove = sum(notMoveLogical(boutStartMLInd:boutEndMLInd));
                numMove = sum(~notMoveLogical(boutStartMLInd:boutEndMLInd));

                % if greater than or equal number of not moving samples,
                %  assign all bouts to not moving
                if (numNotMove >= numMove)
                    notMoveLogical(boutStartMLInd:boutEndMLInd) = true;
                else % otherwise, moving
                    notMoveLogical(boutStartMLInd:boutEndMLInd) = false;
                end
            % one short bout, type change   
            else
                notMoveLogical(boutStartMLInd:boutEndMLInd) = ~boutMoveLogical;
            end
        end
        whichBout = whichBout + 1;
    end
    
    % get start and end indices of not move bouts
    notMoveStartInd = find(diff(notMoveLogical) > 0.9) + 1;
    notMoveEndInd = find(diff(notMoveLogical) < -0.9);

    % add start of trial if fly is not moving at start
    if (notMoveLogical(1))
        notMoveStartInd = [1 notMoveStartInd];
    end
    % add end of trial if fly is not moving at end
    if (notMoveLogical(end))
        notMoveEndInd = [notMoveEndInd length(notMoveLogical)];
    end

    % convert notMoveLogical to indices
    notMoveInd = find(notMoveLogical);

    % make column vectors
    if (isrow(notMoveInd))
        notMoveInd = notMoveInd';
    end
    if (isrow(notMoveStartInd))
        notMoveStartInd = notMoveStartInd';
    end
    if (isrow(notMoveEndInd))
        notMoveEndInd = notMoveEndInd';
    end
end