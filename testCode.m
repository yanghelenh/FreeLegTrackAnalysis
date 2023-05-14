figure; plot(legXShiftRot(40000:41500,refPts.headPtInd),legYShiftRot(40000:41500,refPts.headPtInd), 'x');
hold on;
plot(legXShiftRot(40000:41500,refPts.abdPtInd),legYShiftRot(40000:41500,refPts.abdPtInd), 'x');
plot(legXShiftRot(40000:41500,refPts.thxMidInd),legYShiftRot(40000:41500,refPts.thxMidInd), 'x');


figure; plot(legXShiftRot2(:,refPts.headPtInd),legYShiftRot2(:,refPts.headPtInd), 'x');
hold on;
plot(legXShiftRot2(:,refPts.abdPtInd),legYShiftRot2(:,refPts.abdPtInd), 'x');
plot(legXShiftRot2(:,refPts.thxMidInd),legYShiftRot2(:,refPts.thxMidInd), 'x');

figure; plot(legXShiftRot2(:,refPts.headPtInd),legYShiftRot2(:,refPts.headPtInd), 'x');
hold on;
plot(legXShiftRot2(:,refPts.abdPtInd),legYShiftRot2(:,refPts.abdPtInd), 'x');
plot(legXShiftRot2(:,refPts.thxMidInd),legYShiftRot2(:,refPts.thxMidInd), 'x');
plot(legXShiftRot2(:,1),legYShiftRot2(:,1),'x');


[legTrack.legX, legTrack.legY] = loadTrkFile(trkFilepath);
legY = legTrack.legY(bodytraj.cam.startFr:bodytraj.cam.endFr,:);
legX = legTrack.legX(bodytraj.cam.startFr:bodytraj.cam.endFr,:);
cncX = bodytraj.cnc.x;
cncY = bodytraj.cnc.y;


legY = legTrack.legY;
legX = legTrack.legX;

figure; plot(bodytraj.tZeroed,bodytraj.wrAngSmoM);
% hold on; plot(bodytraj.tZeroed,wrapTo360(flyRawAng + 180));
hold on; 


figure; plot(bodytraj.tZeroed,bodytraj.wrAngSmoM);
% hold on; plot(bodytraj.tZeroed,wrapTo360(flyRawAng + 180));
a = wrapTo360(flyRawAng-90);
a=a(1:length(corrT));
hold on; plot(corrT,a);

plot(t,wrapTo360(flyRawAngAll-90));

flyRawAngAll = flyRawAng;

timeDiff = diff(trialDat.cam.t);
skipInd = find(timeDiff > (mean(timeDiff)) * 3.5);

corrT = trialDat.cam.t;

for i = 1:length(skipInd)
    frameLen = trialDat.cam.t(skipInd(i)+1) - trialDat.cam.t(skipInd(i));

    numFrames = ceil(frameLen / mean(timeDiff));

    addTPts = ((1:(numFrames-1)) * mean(timeDiff)) + trialDat.cam.t(skipInd(i));
    addTPts = addTPts';

    corrT = [corrT(1:skipInd(i)); addTPts; ...
        corrT((skipInd(i)+1):end)];
end

corrT = corrT - corrT(1);



figure; plot(bodytraj.x, bodytraj.y);
hold on;
plot(aLegX(:,10),aLegY(:,10));

figure; plot(bodytraj.tZeroed, bodytraj.x);
hold on;
figure; plot(bodytraj.tZeroed,aLegX(:,1));


testFr = 40000:41500;

for i = 1:6
figure;
plot(legTrack.t(testFr), legTrack.srnLegX(testFr,i));
hold on;
plot(legTrack.t(testFr), legTrack.srnfLegX(testFr,i));

end

figure;
plot(legTrack.t(testFr), legTrack.srnfLegX(testFr,1));
hold on;
plot(legTrack.t(testFr), srnLegSmo(testFr,1));


figure;
plot(legTrack.t(testFr), legVel(testFr,1));

figure;
plot(legTrack.t(testFr), legTrack.aLegX(testFr,1));

for i = 1:13
aLegVelX(:,i) = gradient(legTrack.aLegX(:,i));
end

aLegXVel = findLegVel(legTrack.aLegX, smoParams);

figure;
plot(legTrack.t(testFr), aLegXVel(testFr,1));





whichLeg = 2;

[p, l] = findpeaks(legTrack.srnfLegX(:,whichLeg), 'MinPeakProminence',0.05,...
    'MinPeakDistance',6);

% testFr = 40000:41500;
testFr = 3500:6500;

t = legTrack.t(l);

valInd = find(l>=testFr(1) & l<=testFr(end));
thisT = t(valInd);
thisP = p(valInd);

[p2, l2] = findpeaks(-1*legTrack.srnfLegX(:,whichLeg), 'MinPeakProminence',0.05,...
    'MinPeakDistance',6);



t2 = legTrack.t(l2);

valInd = find(l2>=testFr(1) & l2<=testFr(end));
thisT2 = t2(valInd);
thisP2 = p2(valInd);

figure;
plot(legTrack.t(testFr), legTrack.srnfLegX(testFr,whichLeg));
hold on;
plot(thisT, thisP, 'o');
plot(thisT2, -1*thisP2, 'o');
