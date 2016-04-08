%
% getTempSelective.m
%
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function tempSelective= getTempSelective(nData)
    boxCarWindowLength= 200; % ms
    boxCarWindow      = ones(1,boxCarWindowLength)/(boxCarWindowLength/1000);
    yesActMat         = nData.unit_yes_trial;
    noActMat          = nData.unit_no_trial;
    
    for nAct          = 1:size(yesActMat, 1)
        yesActMat(nAct,:) = conv(yesActMat(nAct,:), boxCarWindow, 'same');
    end
    
    for nAct          = 1:size(noActMat, 1)
        noActMat(nAct,:)  = conv(noActMat(nAct,:), boxCarWindow, 'same');
    end
    
    yesActMat         = yesActMat(:,201:end-200);
    noActMat          = noActMat(:,201:end-200);
    
    [~, tempSelective, ~] = ttest2(yesActMat, noActMat,'dim',1);
    tempSelective(~isnan(tempSelective)) = -log(tempSelective(~isnan(tempSelective)));
    tempSelective(isnan(tempSelective))  = 0;

end