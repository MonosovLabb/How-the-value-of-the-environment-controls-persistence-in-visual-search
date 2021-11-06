% NOTE for users:
%  due to an oddity of Matlab and Excel, when you view the output file
%  in Excel, some of the output columns may have narrow widths and their 
%  numbers may not appear properly (e.g. being rounded to 0 or 1, or being
%  shown as single characters like '#'). However, those columns still 
%  contain the correct numbers. To see them, just resize the width of 
%  those columns in Excel, and the correct numbers should be visible.

clear all;
addpath('functions'); 
addpath('MatlabToolbox-2.8');
% Table calculations
load('Upload_storage_General');
rng(1);
excelOut = 'suppTable.xlsx';
NumberofPerm = 10000;
%% Search Times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars SFRTr SFRTnr mSFRTr mSFRTnr stdSFRTr stdSFRTnr seSFRTr seSFRTnr srocSF spvalueSF seROC_SFRT ...
    RTr RTnr mRTr mRTnr stdRTr stdRTnr seRTr seRTnr srocS spvalueS seROC_RT ...
    rocDiff rocDataIn seROC_SearchDiff

for monkCond = 1:3
    depthPos = monkCond;
    switch monkCond
        case 1
            trialCond = ones(size(storage.monkID)) & storage.goodtrials & ~storage.novelTrial;
        case 2
            trialCond = (storage.monkID == 1) & storage.goodtrials & ~storage.novelTrial;
        case 3
            trialCond = (storage.monkID == 2) & storage.goodtrials & ~storage.novelTrial;
    end  
    tic
    % Single fractal
    SFRTr{depthPos}   = storage.searchDuration(storage.finType==0.3 & storage.rewDur>0 & trialCond, 1);
    SFRTnr{depthPos}  = storage.searchDuration(storage.finType==0.3 & storage.rewDur==0 & trialCond, 1);
    mSFRTr{depthPos}  = mean(SFRTr{depthPos});
    mSFRTnr{depthPos} = mean(SFRTnr{depthPos});
    mdSFRTr{depthPos}  = median(SFRTr{depthPos});
    mdSFRTnr{depthPos} = median(SFRTnr{depthPos});
    stdSFRTr{depthPos} = std(SFRTr{depthPos});
    stdSFRTnr{depthPos} = std(SFRTnr{depthPos});
    seSFRTr{depthPos} = stdSFRTr{depthPos}/length(SFRTr{depthPos});
    seSFRTnr{depthPos} = stdSFRTnr{depthPos}/length(SFRTnr{depthPos});
    [srocSF{depthPos},spvalueSF{depthPos}]=rocarea3(SFRTr{depthPos},SFRTnr{depthPos});
    [CI,BOOTSTAT] = bootci(20000,{@(z) rocarea3(z(z(:,2)==1,1),z(z(:,2)==2,1)),[[SFRTr{depthPos} ; SFRTnr{depthPos}] [1*ones(numel(SFRTr{depthPos}),1); 2*ones(numel(SFRTnr{depthPos}),1)]]},'type','per');
    seROC_SFRT{depthPos} = std(BOOTSTAT);
    toc
    tic
    % Search
    RTr{depthPos}   = storage.searchDuration((storage.finType==0.2 | storage.finType==0.0) & storage.rewDur>0 & trialCond, 1);
    RTnr{depthPos}  = storage.searchDuration((storage.finType==0.2 | storage.finType==0.0) & storage.rewDur==0 & trialCond, 1);
    mRTr{depthPos}  = mean(RTr{depthPos});
    mRTnr{depthPos} = mean(RTnr{depthPos});
    mdRTr{depthPos}  = median(RTr{depthPos});
    mdRTnr{depthPos} = median(RTnr{depthPos});
    stdRTr{depthPos} = std(RTr{depthPos});
    stdRTnr{depthPos} = std(RTnr{depthPos});
    seRTr{depthPos} = stdRTr{depthPos}/length(RTr{depthPos});
    seRTnr{depthPos} = stdRTnr{depthPos}/length(RTnr{depthPos});
    [srocS{depthPos},spvalueS{depthPos}]=rocarea3(RTr{depthPos},RTnr{depthPos});
    [CI,BOOTSTAT] = bootci(20000,{@(z) rocarea3(z(z(:,2)==1,1),z(z(:,2)==2,1)),[[RTr{depthPos} ; RTnr{depthPos}] [1*ones(numel(RTr{depthPos}),1); 2*ones(numel(RTnr{depthPos}),1)]]},'type','per');
    seROC_RT{depthPos} = std(BOOTSTAT);
    toc
    tic
    % Bootstrap of difference between SF and Search
    % Note: rocarea3_area is a wrapper to only return the first value from
    % rocarea, defined at the end of the file
    rocDiff = @(z) rocarea3_area(z(z(:,2)==1,1),z(z(:,2)==2,1)) - rocarea3_area(z(z(:,2)==3,1),z(z(:,2)==4,1));
    rocDataIn = [[SFRTr{depthPos} ; SFRTnr{depthPos}; RTr{depthPos} ; RTnr{depthPos}] [1*ones(numel(SFRTr{depthPos}),1); 2*ones(numel(SFRTnr{depthPos}),1); 3*ones(numel(RTr{depthPos}),1); 4*ones(numel(RTnr{depthPos}),1)]];
    [CIDiff{depthPos},BOOTSTAT] = bootci(20000,{rocDiff, rocDataIn},'alpha',0.001,'type','per');
    if(CI(1) <= 0 && CI(2) >= 0)
        % range does contain 0
        RocDiffSig = 0;
    else
        % range doesn't contain 0, indicate ROC diff is significant
        RocDiffSig = 1;
    end
    seROC_SearchDiff{depthPos} = std(BOOTSTAT);
    toc
    
    % difference
    SearchDiff{depthPos} = abs(mRTr{depthPos} - mRTnr{depthPos});
    SFDiff{depthPos} = abs(mSFRTr{depthPos} - mSFRTnr{depthPos});
    meanDiff = @(x) mean(x(x(:,2)==1,1)) - mean(x(x(:,2)==2,1));
    
    diffData = [SFRTr{depthPos}, ones(size(SFRTr{depthPos}));...
                SFRTnr{depthPos}, 2.*ones(size(SFRTnr{depthPos}))];
    [nrMeanDiffCI{depthPos},BOOTSTAT] = bootci(20000,{meanDiff, rocDataIn},'alpha',0.001,'type','per');
    seSFDiff{depthPos} = std(BOOTSTAT);
    
    diffData = [RTr{depthPos}, ones(size(RTr{depthPos}));...
                RTnr{depthPos}, 2.*ones(size(RTnr{depthPos}))];
    [rMeanDiffCI{depthPos},BOOTSTAT] = bootci(20000,{meanDiff, diffData},'alpha',0.001,'type','per');
    seSearchDiff{depthPos} = std(BOOTSTAT);
end

% Write to excel File
startRowOrig = 3;
startCol = 2; %B
moduleWidth = 7;
for monkCond = 1:3
    rowStart = startRowOrig;
        depthPos = monkCond;
        if monkCond==1
            colStart = char((0)*moduleWidth + 64 + startCol);
            labelStart = char((0)*moduleWidth + 64 + startCol - 1);
        elseif monkCond==2
            colStart = char((2)*moduleWidth + 64 + startCol);
            labelStart = char((2)*moduleWidth + 64 + startCol - 1);
        elseif monkCond==3
            colStart = char((1)*moduleWidth + 64 + startCol);
            labelStart = char((1)*moduleWidth + 64 + startCol - 1);
        end
        
        % monk label
        if monkCond==1
            labelPos = [labelStart int2str(rowStart-2)];
            writematrix("All Data",excelOut,'Sheet',1,'Range',labelPos);
        elseif monkCond==2
            labelPos = [labelStart int2str(rowStart-2)];
            writematrix("Animal M",excelOut,'Sheet',1,'Range',labelPos);
        elseif monkCond==3
            labelPos = [labelStart int2str(rowStart-2)];
            writematrix("Animal B",excelOut,'Sheet',1,'Range',labelPos);
        end
   
        % sf label
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Search duration (s) (Single Fractal Task)",excelOut,'Sheet',1,'Range',labelPos);
        
        % sf meam
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Mean",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([mSFRTr{depthPos};mSFRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % sf stdst
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("SD",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([stdSFRTr{depthPos};stdSFRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % sf sest
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([seSFRTr{depthPos};seSFRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % sf mdst
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Median",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([mdSFRTr{depthPos};mdSFRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % sf n
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("n",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([length(SFRTr{depthPos});length(SFRTnr{depthPos})],excelOut,'Sheet',1,'Range',startPos);
        
        % sf ROC difference RvsNR
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("ROC area (Rew vs No Rew)",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(srocSF{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        rowStart = rowStart;
        labelPos = [labelStart+4 int2str(rowStart)];
        writematrix("p",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart+4 int2str(rowStart)];
        writematrix(spvalueSF{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        rowStart = rowStart;
        labelPos = [labelStart+2 int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart+2 int2str(rowStart)];
        writematrix(seROC_SFRT{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        % Reward vs No Reward Difference
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Difference between Reward and No Reward Mean ",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(SFDiff{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        rowStart = rowStart;
        labelPos = [labelStart+2 int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart+2 int2str(rowStart)];
        writematrix(seSFDiff{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        
        % Spacer row
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix(" ",excelOut,'Sheet',1,'Range',labelPos);
        
        % Search label
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Search duration (s) (Search)",excelOut,'Sheet',1,'Range',labelPos);
        
        % search m
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Mean",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([mRTr{depthPos};mRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % search stdst
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("SD",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([stdRTr{depthPos};stdRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % search sest
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([seRTr{depthPos};seRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % search mdst
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Median",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([mdRTr{depthPos};mdRTnr{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % search n
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("n",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([length(RTr{depthPos});length(RTnr{depthPos})],excelOut,'Sheet',1,'Range',startPos);
        
        % search ROC difference RvsNR
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("ROC area (Rew vs No Rew",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(srocS{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        rowStart = rowStart;
        labelPos = [labelStart+4 int2str(rowStart)];
        writematrix("p",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart+4 int2str(rowStart)];
        writematrix(spvalueS{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        rowStart = rowStart;
        labelPos = [labelStart+2 int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart+2 int2str(rowStart)];
        writematrix(seROC_RT{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        % Reward vs No Reward Difference
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Difference between Reward and No Reward Mean ",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(SearchDiff{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        rowStart = rowStart;
        labelPos = [labelStart+2 int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart+2 int2str(rowStart)];
        writematrix(seSearchDiff{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        % Spacer row
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix(" ",excelOut,'Sheet',1,'Range',labelPos);
        
        % % Bootstrap of difference in ROC area between SF and Search -
        % diff
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("ROC area difference across tasks",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix([srocS{depthPos} - srocSF{depthPos}],excelOut,'Sheet',1,'Range',startPos);
        
        % upper and lower bound - excludes zero test
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Bootstrap 99.9% CI",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(CIDiff{depthPos}*-1,excelOut,'Sheet',1,'Range',startPos); % -1 multiplier to get into the same format as the main difference...
        
        % Spacer row
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix(" ",excelOut,'Sheet',1,'Range',labelPos);
end

%% Search Times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except storage excelOut NumberofPerm rowStart
SRT_calc_NR = {};
SRT_calc_R = {};

for monkCond = 1:3
    for envCond = 1:3
        depthPos = (monkCond-1)*3 + envCond;
        switch monkCond
            case 1
                trialCond = ones(size(storage.monkID)) & storage.goodtrials & ~storage.novelTrial;
            case 2
                trialCond = (storage.monkID == 1) & storage.goodtrials & ~storage.novelTrial;
            case 3
                trialCond = (storage.monkID == 2) & storage.goodtrials & ~storage.novelTrial;
        end
        switch envCond
            case 1
                % no restrictions
            case 2
                trialCond = trialCond & storage.envID < 6;
            case 3
                trialCond = trialCond & storage.envID > 6;
        end
        SRT_calc_NR{depthPos} = [];
        SRT_calc_R{depthPos} = [];
        for i = 1:5
            env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & trialCond;
            % NO Rew
            tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'==0);
            st{i,1,depthPos} =  storage.searchDuration(tempIndex,1);
            SRT_calc_NR{depthPos} = [SRT_calc_NR{depthPos}; [st{i,1,depthPos}, repmat(i, size(st{i,1,depthPos},1), 1)]];
            % Rew
            tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'>0);
            st{i,2,depthPos} =  storage.searchDuration(tempIndex,1);
            SRT_calc_R{depthPos} =  [SRT_calc_R{depthPos}; [st{i,2,depthPos}, repmat(i, size(st{i,2,depthPos},1), 1)]];
        end        
    end
end

% ST
% MeasID: FSOH39
% MeasID: MVV69A 
mst = cellfun(@mean,st);
stdst = cellfun(@std,st);
sest = stdst./sqrt(cellfun(@length,st));
mdst = cellfun(@median,st);
% n
nst = cellfun(@length,st);

% Significance analyses
for monkCond = 1:3
    for envCond = 1:3
        depthPos = (monkCond-1)*3 + envCond;
        % Spearman Correlation
        [pvalNRrank, rNRrank] = permutation_pair_test_fast(SRT_calc_NR{depthPos}(:,1),SRT_calc_NR{depthPos}(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
        [pvalRrank,  rRrank ] = permutation_pair_test_fast(SRT_calc_R{depthPos}(:,1),SRT_calc_R{depthPos}(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
        tempSpearman{depthPos} = [pvalRrank,  rRrank; pvalNRrank, rNRrank];
    end
end

nst(1,2,:) = nan;
nst(5,1,:) = nan;

for cleanIter = 1:9
    pstor{cleanIter}(1,2) = nan;
    pstor{cleanIter}(4,1) = nan;
end

% Write to excel File
startRowOrig = rowStart+4;
startCol = 2; %B
moduleWidth = 7;
for monkCond = 1:3
    rowStart = startRowOrig;
    for envCond = 1:3
        depthPos = (monkCond-1)*3 + envCond;
        if monkCond==1
            colStart = char((0)*moduleWidth + 64 + startCol);
            labelStart = char((0)*moduleWidth + 64 + startCol - 1);
        elseif monkCond==2
            colStart = char((2)*moduleWidth + 64 + startCol);
            labelStart = char((2)*moduleWidth + 64 + startCol - 1);
        elseif monkCond==3
            colStart = char((1)*moduleWidth + 64 + startCol);
            labelStart = char((1)*moduleWidth + 64 + startCol - 1);
        end
        
        labelPos = [labelStart int2str(rowStart)];
        if envCond==1
            writematrix("Search duration (s) (Search Task)",excelOut,'Sheet',1,'Range',labelPos);
        elseif envCond==2
            writematrix("Search duration (s) (Search Task, Environment Set 2 only)",excelOut,'Sheet',1,'Range',labelPos);
        else
            writematrix("Search duration (s) (Search Task, Environment Set 1 only)",excelOut,'Sheet',1,'Range',labelPos);
        end

        % mst
        rowStart = rowStart+1;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Mean",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(flipud(mst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
        
        % stdst
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("SD",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(flipud(stdst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
        
        % sest
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(flipud(sest(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
        
        % mdst
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Median",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(flipud(mdst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
        
        % n
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("n",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(flipud(nst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
                
        % Spearman's (rankcorr)
        rowStart = rowStart+2;
        labelPos = [labelStart int2str(rowStart)];
        writematrix("Corr.w/P(rew)",excelOut,'Sheet',1,'Range',labelPos);
        startPos = [colStart int2str(rowStart)];
        writematrix(tempSpearman{depthPos},excelOut,'Sheet',1,'Range',startPos);
        
        % Spacer row
        rowStart = rowStart+4;
        labelPos = [labelStart int2str(rowStart)];
        writematrix(" ",excelOut,'Sheet',1,'Range',labelPos);
    end
end

%% Persistent Search Duration
SRT_calc_NR = {};
SRT_calc_R = {};
st = cell(5,2,3);
for monkCond = 1:3
    depthPos = monkCond;
    switch monkCond
        case 1
            trialCond = ones(size(storage.monkID));
        case 2
            trialCond = storage.monkID == 1;
        case 3
            trialCond = storage.monkID == 2;
    end
    SRT_calc_NR{depthPos} = [];
    SRT_calc_R{depthPos} = [];

    for i = 1:5
        env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & trialCond;
        % NO Rew
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'==0 & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
        st{i,1,depthPos} =  storage.firstSacc2TargPersSearchDur(tempIndex)';
        st{i,1,depthPos} =  st{i,1,depthPos}(~isnan(st{i,1,depthPos}));
        SRT_calc_NR{depthPos} = [SRT_calc_NR{depthPos}; [st{i,1,depthPos}, repmat(i, size(st{i,1,depthPos},1), 1)]];
        % Rew
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'>0 & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
        st{i,2,depthPos} =  storage.firstSacc2TargPersSearchDur(tempIndex)';
        st{i,2,depthPos} =  st{i,2,depthPos}(~isnan(st{i,2,depthPos}));
        SRT_calc_R{depthPos} =  [SRT_calc_R{depthPos}; [st{i,2,depthPos}, repmat(i, size(st{i,2,depthPos},1), 1)]];
    end
end

% ST
mst = cellfun(@mean,st);
stdst = cellfun(@std,st);
sest = stdst./sqrt(cellfun(@length,st));
mdst = cellfun(@median,st);

% Significance analyses
for monkCond = 1:3
    depthPos = monkCond;
    % Spearman Correlation
    [pvalNRrank, rNRrank] = permutation_pair_test_fast(SRT_calc_NR{depthPos}(:,1),SRT_calc_NR{depthPos}(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
    [pvalRrank,  rRrank ] = permutation_pair_test_fast(SRT_calc_R{depthPos}(:,1),SRT_calc_R{depthPos}(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
    tempSpearman{depthPos} = [pvalRrank,  rRrank; pvalNRrank, rNRrank];
end

for cleanIter = 1:9
    pstor{cleanIter}(1,2) = nan;
    pstor{cleanIter}(4,1) = nan;
end

% Write to excel File
startRowOrig = rowStart+4;
startCol = 2; %B
moduleWidth = 7;
for monkCond = 1:3
    rowStart = startRowOrig;
    depthPos = monkCond;
    if monkCond==1
        colStart = char((0)*moduleWidth + 64 + startCol);
        labelStart = char((0)*moduleWidth + 64 + startCol - 1);
    elseif monkCond==2
        colStart = char((2)*moduleWidth + 64 + startCol);
        labelStart = char((2)*moduleWidth + 64 + startCol - 1);
    elseif monkCond==3
        colStart = char((1)*moduleWidth + 64 + startCol);
        labelStart = char((1)*moduleWidth + 64 + startCol - 1);
    end
    
    labelPos = [labelStart int2str(rowStart)];
    writematrix("Remaining search duration after gazing at target with a direct saccade (s) (Search Task)",excelOut,'Sheet',1,'Range',labelPos);
    
    % mst
    rowStart = rowStart+1;
    labelPos = [labelStart int2str(rowStart)];
    writematrix("Mean",excelOut,'Sheet',1,'Range',labelPos);
    startPos = [colStart int2str(rowStart)];
    writematrix(flipud(mst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
    
    % stdst
    rowStart = rowStart+2;
    labelPos = [labelStart int2str(rowStart)];
    writematrix("SD",excelOut,'Sheet',1,'Range',labelPos);
    startPos = [colStart int2str(rowStart)];
    writematrix(flipud(stdst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
    
    % sest
    rowStart = rowStart+2;
    labelPos = [labelStart int2str(rowStart)];
    writematrix("SE",excelOut,'Sheet',1,'Range',labelPos);
    startPos = [colStart int2str(rowStart)];
    writematrix(flipud(sest(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);
    
    % mdst
    rowStart = rowStart+2;
    labelPos = [labelStart int2str(rowStart)];
    writematrix("Median",excelOut,'Sheet',1,'Range',labelPos);
    startPos = [colStart int2str(rowStart)];
    writematrix(flipud(mdst(:,:,depthPos)'),excelOut,'Sheet',1,'Range',startPos);

    % Spearman's (rankcorr)
    rowStart = rowStart+2;
    labelPos = [labelStart int2str(rowStart)];
    writematrix("Corr.w/P(rew)",excelOut,'Sheet',1,'Range',labelPos);
    startPos = [colStart int2str(rowStart)];
    writematrix(tempSpearman{depthPos},excelOut,'Sheet',1,'Range',startPos);
    
    % Spacer row
    rowStart = rowStart+4;
    labelPos = [labelStart int2str(rowStart)];
    writematrix(" ",excelOut,'Sheet',1,'Range',labelPos);
end

%% ERRORS Search timeout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except storage NumberofPerm excelOut rowStart

for i = 1:5
    env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & ~storage.novelTrial & storage.goodtrials;

    % NO Rew
    n{i,1} = (logical(env(:,i) & (storage.finType'==3.0 | storage.finType'==0.2) & storage.rewDur'==0));
    SearchTimeout{i,1} = (logical(env(:,i) & storage.finType'==3.0 & storage.rewDur'==0));
    
    % Rew
    n{i,2} = (logical(env(:,i) & (storage.finType'==3.0 | storage.finType'==0.2) & storage.rewDur'>0));
    SearchTimeout{i,2} = (logical(env(:,i) & storage.finType'==3.0 & storage.rewDur'>0));
end

SRT_calc_NR = {};
SRT_calc_R = {};
st = cell(5,2,3);
for monkCond = 1:3
    depthPos = monkCond;
    switch monkCond
        case 1
            trialCond = ones(size(storage.monkID));
        case 2
            trialCond = storage.monkID == 1;
        case 3
            trialCond = storage.monkID == 2;
    end
    SRT_calc_NR{depthPos} = [];
    SRT_calc_R{depthPos} = [];

    for i = 1:4
        % NO Rew
        st{i,1,depthPos} =  sum(SearchTimeout{i,1} & trialCond')/sum(n{i,1});
        SRT_calc_NR{depthPos} = [SRT_calc_NR{depthPos}; [st{i,1,depthPos}, repmat(i, size(st{i,1,depthPos},1), 1)]];
        % Rew
        st{i+1,2,depthPos} =  sum(SearchTimeout{i+1,2} & trialCond')/sum(n{i+1,2});
        SRT_calc_R{depthPos} =  [SRT_calc_R{depthPos}; [st{i+1,2,depthPos}, repmat(i+1, size(st{i+1,2,depthPos},1), 1)]];
    end
end

clear mst
for iter = 1:4
    mst{1}(iter,1) = st{iter,1,1};
    mst{1}(iter+1,2) = st{iter+1,2,1};
    mst{2}(iter,1) = st{iter,1,2};
    mst{2}(iter+1,2) = st{iter+1,2,2};
    mst{3}(iter,1) = st{iter,1,3};
    mst{3}(iter+1,2) = st{iter+1,2,3};
end
for iter = 1:3
    mst{iter} = mst{iter}';
    mst{iter}(2,1) = nan;
    mst{iter}(1,5) = nan;
end

% Write to excel File
startRowOrig = rowStart+9;
startCol = 2; %B
moduleWidth = 7;
for monkCond = 1:3
    rowStart = startRowOrig;
    depthPos = monkCond;
    if monkCond==1
        colStart = char((0)*moduleWidth + 64 + startCol);
        labelStart = char((0)*moduleWidth + 64 + startCol - 1);
    elseif monkCond==2
        colStart = char((2)*moduleWidth + 64 + startCol);
        labelStart = char((2)*moduleWidth + 64 + startCol - 1);
    elseif monkCond==3
        colStart = char((1)*moduleWidth + 64 + startCol);
        labelStart = char((1)*moduleWidth + 64 + startCol - 1);
    end
    
    labelPos = [labelStart int2str(rowStart)];
    writematrix("Probability of search errors (Search Task)",excelOut,'Sheet',1,'Range',labelPos);
    
    % mst
    rowStart = rowStart+1;
    labelPos = [labelStart int2str(rowStart)];
    writematrix("Mean",excelOut,'Sheet',1,'Range',labelPos);
    startPos = [colStart int2str(rowStart)];
    writematrix(flipud(mst{depthPos}),excelOut,'Sheet',1,'Range',startPos);
end




function retVal = rocarea3_area(a, b) % return only the first value from rocarea3 for area difference calculations
    temp = rocarea3(a,b,'exact');
    retVal = temp(1);
end



function retVal = corrMod(a) % return only the first value from rocarea3 for area difference calculations
    temp = corr(a, 'type', 'spearman');
    retVal = temp(1,2);
end
