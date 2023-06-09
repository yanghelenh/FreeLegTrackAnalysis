% 6/8/23 - Noticed that heading could be off by 180 deg, which made the fly
%  appear to walk backwards when it wasn't. Manually went through all 2018
%  trials (doesn't affect 2022, but those have other issues), and corrected
%  it using this script.

figure;
plot(bodytraj.tZeroed, bodytraj.unWrAng);

%%
figure; plot(legTrack.t, legTrack.flyRawAng); 
hold on; 
plot(bodytraj.tZeroed, wrapTo360(bodytraj.unWrAng));

%%
figure; plot(bodytraj.tZeroed, bodytraj.fwdVelSmoS);

%%
figure; plot(legTrack.t, legTrack.flyRawAng); 
hold on; 
plot(bodytraj.tZeroed, wrapTo360(bodytraj.unWrAng+180));

%%
bodytraj.unWrAng(1:3659) = bodytraj.unWrAng(1:3659) + 180;

%%

bodytraj.unWrAng = bodytraj.unWrAng + 180;

            %%

                            % unwrap angles        
                jumpThresh = 120;
                unWrAngs = [];
                angOffset = 0;
        
                for j = 1:length(bodytraj.rawAngs)
                    curAng = bodytraj.rawAngs(j);
                    if (j == 1)
                        unWrAngs = [unWrAngs; curAng];
                        continue;
                    end
        
                    curAng = curAng + angOffset;
        
                    if ((unWrAngs(end) - curAng) > jumpThresh)
                        angOffset = angOffset + 180;
                        curAng = curAng + 180;
                    elseif ((curAng - unWrAngs(end)) > jumpThresh)
                        angOffset = angOffset - 180;
                        curAng = curAng - 180;
                    end
        
                    unWrAngs = [unWrAngs; curAng];
                end
                
                % correct angles from ellipse fitting definition, camera flip
                unWrAngs = 90 - unWrAngs;
                unWrAngs = -1 * unWrAngs;

                bodytraj.unWrAng = unWrAngs;

%%

    % define constants
    % max number of adjacent fly not present frames to smooth over
    MAX_FLYPRESENT_FRAMES = 13;
    % smoothing of fly behavior, window for smoothdata function, 
    %  method gaussian
    SIGMA_MILD = 15; 
    SIGMA_STRONG = 50; % aggressive smoothing

           % smoothing on position
            if (SIGMA_MILD ~= 0)
                % mild smoothing
                bodytraj.xSmoM = smoothdata(bodytraj.x, 'gaussian', ...
                    SIGMA_MILD);
                bodytraj.ySmoM = smoothdata(bodytraj.y, 'gaussian', ...
                    SIGMA_MILD);
                bodytraj.unWrAngSmoM = smoothdata(bodytraj.unWrAng, ...
                    'gaussian', SIGMA_MILD);
                % aggressive smothing
                bodytraj.xSmoS = smoothdata(bodytraj.x, 'gaussian', ...
                    SIGMA_STRONG);
                bodytraj.ySmoS = smoothdata(bodytraj.y, 'gaussian', ...
                    SIGMA_STRONG);
                bodytraj.unWrAngSmoS = smoothdata(bodytraj.unWrAng, ...
                    'gaussian', SIGMA_STRONG);
            end
            
            % get wrapped angles
            bodytraj.wrAng = wrapTo360(bodytraj.unWrAng); 
            bodytraj.wrAngSmoM = wrapTo360(bodytraj.unWrAngSmoM);
            bodytraj.wrAngSmoS = wrapTo360(bodytraj.unWrAngSmoS);
            
            % compute velocities, on smoothed data
            % mild smoothing
            xVelM = gradient(bodytraj.xSmoM, bodytraj.tZeroed);
            yVelM = gradient(bodytraj.ySmoM, bodytraj.tZeroed);
            % fly's angular velocity, turns to fly's own right are positive
            angVelM = gradient(bodytraj.unWrAngSmoM, bodytraj.tZeroed);

            % strong smoothing
            xVelS = gradient(bodytraj.xSmoS, bodytraj.tZeroed);
            yVelS = gradient(bodytraj.ySmoS, bodytraj.tZeroed);
            % fly's angular velocity, turns to fly's own right are positive
            angVelS = gradient(bodytraj.unWrAngSmoS, bodytraj.tZeroed);
            
            % smooth velocities
            bodytraj.xVelSmoM = smoothdata(xVelM, 'gaussian', SIGMA_MILD);
            bodytraj.yVelSmoM = smoothdata(yVelM, 'gaussian', SIGMA_MILD);
            bodytraj.angVelSmoM = smoothdata(angVelM, 'gaussian', SIGMA_MILD);

            bodytraj.xVelSmoS = smoothdata(xVelS, 'gaussian', SIGMA_STRONG);
            bodytraj.yVelSmoS = smoothdata(yVelS, 'gaussian', SIGMA_STRONG);
            bodytraj.angVelSmoS = smoothdata(angVelS, 'gaussian',...
                SIGMA_STRONG);


            % translational velocity in any direction
            bodytraj.transVelSmoM = sqrt(bodytraj.xVelSmoM.^2 + ...
                bodytraj.yVelSmoM.^2);
            bodytraj.transVelSmoS = sqrt(bodytraj.xVelSmoS.^2 + ...
                bodytraj.yVelSmoS.^2);
            
            % forward velocity, mild smoothing
            velAngM = atan2d(bodytraj.yVelSmoM, bodytraj.xVelSmoM);
            angDiffM = velAngM - bodytraj.wrAngSmoM; 
            % fly's forward velocity, forward is positive
            bodytraj.fwdVelSmoM = bodytraj.transVelSmoM .* cosd(angDiffM);

            % forward velocity, strong smoothing
            velAngS = atan2d(bodytraj.yVelSmoS, bodytraj.xVelSmoS);
            angDiffS = velAngS - bodytraj.wrAngSmoS; 
            % fly's forward velocity, forward is positive
            bodytraj.fwdVelSmoS = bodytraj.transVelSmoS .* cosd(angDiffS);
            
            % fly's lateral velocity, movement to fly's own right is positive
            % mild smoothing
            bodytraj.latVelSmoM = bodytraj.transVelSmoM .* sind(angDiffM); 
            % strong smoothing
            bodytraj.latVelSmoS = bodytraj.transVelSmoS .* sind(angDiffS);  
            
            % distance covered (only mild smoothing)
            xCovered = sum(abs(diff(bodytraj.xSmoM)));
            yCovered = sum(abs(diff(bodytraj.ySmoM)));
            bodytraj.distCovered = xCovered + yCovered;
            bodytraj.angCovered = sum(abs(diff(bodytraj.unWrAngSmoM)));
            
            % duration of trial
            bodytraj.duration = bodytraj.t(end) - bodytraj.t(1);




%%

thisPDataName = 'trial-31-20181108-182724.mat';
save(thisPDataName, 'bodytraj', '-append');





