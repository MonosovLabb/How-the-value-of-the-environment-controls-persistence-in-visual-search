clear all; %close all; clc;
addpath('functions');
addpath('MatlabToolbox-2.8');
load('Upload_storage_General.mat');


rng(1);

excelOut = 'tableTemp.xlsx';

nrewColor = [0.6235    0.6196    0.9961];
rewColor  = [0.9686    0.4627    0.5137];
meanColor = [1.0000    0.9961    0.3294];
sizeMeanMark = 20;
astFontsize = 24;
textFontsize = 12;

NumberofPerm=10000; %shuffle 10,000 times to get P-values for correlations

%% Panel 1: Single Fractal, Search means; ROCs
figure; 

% Single fractal
SFRTr   = storage.searchDuration(storage.finType==0.3 & storage.rewDur>0 & ~storage.novelTrial & storage.goodtrials, 1);
SFRTnr  = storage.searchDuration(storage.finType==0.3 & storage.rewDur==0 & ~storage.novelTrial & storage.goodtrials, 1);
mSFRTr  = mean(SFRTr);
mSFRTnr = mean(SFRTnr);
stdSFRTr = std(SFRTr);
stdSFRTnr = std(SFRTnr);
seSFRTr = stdSFRTr/length(SFRTr);
seSFRTnr = stdSFRTnr/length(SFRTnr);
[srocSF,spvalueSF]=rocarea3(SFRTr,SFRTnr);
BOOTSTAT = bootstrp(20000, @rocWrapper,[[SFRTr ; SFRTnr] [1*ones(numel(SFRTr),1); 2*ones(numel(SFRTnr),1)]]);
seROC_SFRT = std(BOOTSTAT);

% Search
RTr   = storage.searchDuration((storage.finType==0.2 | storage.finType==0.0) & storage.rewDur>0 & ~storage.novelTrial & storage.goodtrials, 1);
RTnr  = storage.searchDuration((storage.finType==0.2 | storage.finType==0.0) & storage.rewDur==0 & ~storage.novelTrial & storage.goodtrials, 1);
mRTr  = mean(RTr);
mRTnr = mean(RTnr);
stdRTr = std(RTr);
stdRTnr = std(RTnr);
seRTr = stdRTr/length(RTr);
seRTnr = stdRTnr/length(RTnr);
[srocS,spvalueS]=rocarea3(RTr,RTnr);
BOOTSTAT = bootstrp(20000, @rocWrapper,[[RTr ; RTnr] [1*ones(numel(RTr),1); 2*ones(numel(RTnr),1)]]);
seROC_RT = std(BOOTSTAT);

% Bootstrap of difference between SF and Search
% Note: rocarea3_area is a wrapper to only return the first value from
% rocarea, defined at the end of the file
rocDiff = @(z) rocarea3_area(z(z(:,2)==1,1),z(z(:,2)==2,1)) - rocarea3_area(z(z(:,2)==3,1),z(z(:,2)==4,1));
rocDataIn = [[SFRTr ; SFRTnr; RTr ; RTnr] [1*ones(numel(SFRTr),1); 2*ones(numel(SFRTnr),1); 3*ones(numel(RTr),1); 4*ones(numel(RTnr),1)]];
[CI,BOOTSTAT] = bootci(20000,{rocDiff, rocDataIn},'alpha',0.001,'type','per');
if(CI(1) <= 0 && CI(2) >= 0)
    % range does contain 0
    RocDiffSig = 0;
else
    % range doesn't contain 0, indicate ROC diff is significant
    RocDiffSig = 1;
end

% Mean Time Plots
xs = 1:4;
linew = 1;
offset = 0.05;
rewFill = rewColor; 
nRewFill = nrewColor;
iosr.statistics.boxPlot(SFRTnr, 'showviolin', false, 'boxcolor', nRewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', xs(1)+offset, 'groupWidth', 1.2, 'showMean', true, 'meanColor',  'b', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;
iosr.statistics.boxPlot(SFRTr, 'showviolin', false, 'boxcolor', rewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', xs(2)-offset, 'groupWidth', 1.2, 'showMean', true, 'meanColor',  'r', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;
iosr.statistics.boxPlot(RTnr, 'showviolin', false, 'boxcolor', nRewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', xs(3)+offset, 'groupWidth', 1.2, 'showMean', true, 'meanColor',  'b', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;
iosr.statistics.boxPlot(RTr, 'showviolin', false, 'boxcolor', rewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', xs(4)-offset, 'groupWidth', 1.2, 'showMean', true, 'meanColor',  'r', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;
set(gca, 'XTick', [1.5 3.5], 'XTickLabel', {'Single Fractal', 'Search'});

ylabel('Search Duration (s)');
tempn = sum([length(SFRTnr) length(SFRTr) length(RTnr) length(RTr)]);
currY = ylim;
xlim([0.5 4.5]);
ylim([0 2.5])

rvnrSF = ranksum(SFRTr, SFRTnr);

if rvnrSF<0.001
    text(1.5, max(mSFRTr, mSFRTnr)+0.1, '***','FontSize',astFontsize, 'HorizontalAlignment', 'center');
elseif rvnrSF<0.05
    text(1.5, max(mSFRTr, mSFRTnr)+0.1, '*','FontSize',astFontsize, 'HorizontalAlignment',  'center');
else
    text(1.5, max(mSFRTr, mSFRTnr)+0.1, 'n.s.','FontSize',textFontsize, 'HorizontalAlignment',  'center');
end

rvnrS = ranksum(RTr, RTnr);

if rvnrS<0.001
    text(3.5, max(mRTr, mRTnr)+0.1, '***','FontSize',astFontsize, 'HorizontalAlignment',  'center');
elseif rvnrS<0.05
    text(3.5, max(mRTr, mRTnr)+0.1, '*','FontSize',astFontsize, 'HorizontalAlignment',  'center');
else
    text(3.5, max(mRTr, mRTnr)+0.1, 'n.s.','FontSize',textFontsize, 'HorizontalAlignment',  'center');
end


%% Discriminability Plots
figure;

mydata = [srocSF, srocS];
b = bar(1:2, mydata); hold on;
set(b,'FaceColor','w');
colorErr = 'k'; 
errorbar(mydata, [seROC_SFRT, seROC_RT], '.', 'Color', colorErr, 'XData', 1:2);
set(gca, 'XTick', [1:2], 'XTickLabel', {'Single Fractal', 'Search'});

if RocDiffSig
    text(1.5, max(srocSF,srocS)+0.15, '***','FontSize',astFontsize, 'HorizontalAlignment', 'center');
else
    text(1.5, srocSF+0.15, 'n.s.','FontSize',textFontsize, 'HorizontalAlignment',  'center');
end

if spvalueSF<0.00001
    text(1, srocSF+0.1, '***','FontSize',astFontsize, 'HorizontalAlignment', 'center');
elseif spvalueSF<0.05
    text(1, srocSF+0.1, '*','FontSize',astFontsize, 'HorizontalAlignment',  'center');
else
    text(1, srocSF+0.1, 'n.s.','FontSize',textFontsize, 'HorizontalAlignment',  'center');
end

if spvalueS<0.00001
    text(2, srocS+0.1, '***','FontSize',astFontsize, 'HorizontalAlignment',  'center');
elseif spvalueS<0.05
    text(2, srocS+0.1, '*','FontSize',astFontsize, 'HorizontalAlignment',  'center');
else
    text(2, srocS+0.1, 'n.s.','FontSize',textFontsize, 'HorizontalAlignment',  'center');
end

ylim([0.5 1]);
xlim([0.25 2.75])
ref = refline(0,0.5);
ref.Color = rewColor;
set(gca, 'YTick', [0:0.25:1], 'YTickLabel', {'0', '0.25','0.5','0.75','1.00'});
ylabel({'Reward vs no reward discriminability','(ROC area)'});

%% Panel 2 C - Search Means 
figure;

xOffset = 0.1;
noRewXs = 2*(1:4)-1+xOffset;
rewXs = 2*((1:4))+2 - xOffset;

SRT_calc_NR = [];
SRT_calc_R  = [];
for i = 1:5
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    % NO Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'==0 & ~storage.novelTrial' & storage.goodtrials');
    st{i,1,1} =  storage.searchDuration(tempIndex,1);
    SRT_calc_NR = [SRT_calc_NR; [st{i,1,1}, repmat(i, size(st{i,1,1},1), 1)]];
    % Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'>0 & ~storage.novelTrial' & storage.goodtrials');
    st{i,2,1} =  storage.searchDuration(tempIndex,1);
    SRT_calc_R =  [SRT_calc_R; [st{i,2,1}, repmat(i, size(st{i,2,1},1), 1)]];
end

% ST
mst = cellfun(@median,st);
stdst = cellfun(@std,st);
sest = stdst./sqrt(cellfun(@length,st));

% Significance analyses
for i = 1:3
    pstor(i,1) = ranksum(st{i,1,1}, st{i+1,1,1});
    pstor(i+1,2) = ranksum(st{i+1,2,1}, st{i+2,2,1});
end
[pvalNR, rNR] = permutation_pair_test_fast(SRT_calc_NR(:,1),SRT_calc_NR(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
[pvalR,  rR ] = permutation_pair_test_fast(SRT_calc_R(:,1),SRT_calc_R(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear

rewFill = rewColor;
nRewFill = nrewColor;
linew = 1;
% Mean Search time & STD
for overlap = 1%:2
    hold on;
    for i = 1:4
        % NoRew
        iosr.statistics.boxPlot(st{i,1}, 'showviolin', false, 'boxcolor', nRewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', noRewXs(i), 'boxWidth', 0.8, 'showMean', true, 'meanColor',  'b', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;

        % Rew
        iosr.statistics.boxPlot(st{i+1,2}, 'showviolin', false, 'boxcolor', rewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', rewXs(i), 'boxWidth', 0.8, 'showMean', true, 'meanColor',  'r', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;
    end
    
    % Plot p-values on lines
    xlim([0.5 10.5])
    currY = ylim;
    ylim([0 currY(2)])
    xlabel('Environment reward probability')
    ylabel('Search Duration (s)')
    set(gca, 'XTick', [1.5, 3.5, 5.5, 7.5, 9.5, 11.5], 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
    
    clear net sing
    [net, sing] = chartCalcs(st, NumberofPerm);
    textPrintMeans(st, net, sing, noRewXs, rewXs, 0.01, nrewColor, rewColor);
    
    xes = xlim;
    hold off
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
end


%% 2D: Correlations
clearvars -except storage figNum NumberofPerm indiv nrewColor rewColor
figure;

% Search NR
SRTnr = [];
SRTnr_calc = [];
env = [];
for i = 1:4
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    index = env(:,i) & storage.rewDur'==0 & (storage.finType'==0.0 | storage.finType'==0.2) & ~storage.novelTrial' & storage.goodtrials';
    RTs = storage.searchDuration(index, 1);
    SRTnr_calc = [SRTnr_calc; [RTs, repmat(i, length(RTs), 1)]];
end
tic
[pval_S_nr, r_S_nr]=permutation_pair_test_fast(SRTnr_calc(:,1),SRTnr_calc(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear

[R,P,confLow{1},confHigh{1}] = corrcoef(SRTnr_calc(:,1),SRTnr_calc(:,2));
nTemp = length(SRTnr_calc(:,1));

% Search R
SRTr = [];
SRTr_calc = [];
env = [];
for i = 2:5
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    index = env(:,i) & storage.rewDur'>0 & (storage.finType'==0.0 | storage.finType'==0.2) & ~storage.novelTrial' & storage.goodtrials';
    RTs = storage.searchDuration(index, 1);
    SRTr_calc = [SRTr_calc; [RTs, repmat(i, length(RTs), 1)]];
end
[pval_S_r, r_S_r]=permutation_pair_test_fast(SRTr_calc(:,1),SRTr_calc(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
[R,P,confLow{2},confHigh{2}] = corrcoef(SRTr_calc(:,1),SRTr_calc(:,2));
nTemp = nTemp + length(SRTr_calc(:,1));

searchNR = SRTnr_calc;
searchR = SRTr_calc;

% Single fractal STs for those target fractals associated with a single environment
envTargs = [87, 88, 89, 90;
    84, 85, 86, 91;
    82, 83, 92, 93;
    81, 94, 95, 96;
    97, 98, 99, 100];

% Un-precued SF NR
SRTnr = [];
SRTnr_calc = [];
env = [];
for i = 1:4
    env(:,i) = storage.TargID' == envTargs(i,1) | storage.TargID' == envTargs(i,2) | storage.TargID' == envTargs(i,3) | storage.TargID' == envTargs(i,4);
    index = env(:,i) & storage.rewDur'==0 & (storage.finType'==0.3) & ~storage.novelTrial' & storage.goodtrials';
    RTs = storage.searchDuration(index, 1);
    SRTnr_calc = [SRTnr_calc; [RTs, repmat(i, length(RTs), 1)]];
end
zz = SRTnr;
[pval_SF_nr, r_SF_nr]=permutation_pair_test_fast(SRTnr_calc(:,1),SRTnr_calc(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
[R,P,confLow{3},confHigh{3}] = corrcoef(SRTnr_calc(:,1),SRTnr_calc(:,2));
nTemp = nTemp + length(SRTnr_calc(:,1));

% Un-precued SF R
SRTr = [];
SRTr_calc = [];
env = [];
for i = 2:5
    env(:,i) = storage.TargID' == envTargs(i,1) | storage.TargID' == envTargs(i,2) | storage.TargID' == envTargs(i,3) | storage.TargID' == envTargs(i,4);
    index = env(:,i) & storage.rewDur'>0 & (storage.finType'==0.3) & ~storage.novelTrial' & storage.goodtrials';
    RTs = storage.searchDuration(index, 1);
    SRTr_calc = [SRTr_calc; [RTs, repmat(i, length(RTs), 1)]];
end

zz = SRTr;
[pval_SF_r, r_SF_r]=permutation_pair_test_fast(SRTr_calc(:,1),SRTr_calc(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
[R,P,confLow{4},confHigh{4}] = corrcoef(SRTr_calc(:,1),SRTr_calc(:,2));
nTemp = nTemp + length(SRTr_calc(:,1));

singfracNR = SRTnr_calc;
singfracR = SRTr_calc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wipe Current Storage & Import Precued %%%%%%%%%%%%%%%%%%%%%%%%%%
clear storage;
load('Upload_storage_Fig2Precued.mat')

% Precued SF NR
SRTnr = [];
SRTnr_calc = [];
env = [];
for i = 1:4
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    index = env(:,i) & storage.rewDur'==0 & storage.finType'==0.3 & storage.goodtrials';
    RTs = storage.searchDuration(index, 1);
    SRTnr_calc = [SRTnr_calc; [RTs, repmat(i, length(RTs), 1)]];
end
[pval_SFpq_nr, r_SFpq_nr]=permutation_pair_test_fast(SRTnr_calc(:,1),SRTnr_calc(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
nTemp = nTemp + length(SRTnr_calc(:,1));

% Precued SF R
SRTr = [];
SRTr_calc = [];
env = [];
for i = 2:5
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    index = env(:,i) & storage.rewDur'>0 & storage.finType'==0.3 & storage.goodtrials';
    RTs = storage.searchDuration(index, 1);
    SRTr_calc = [SRTr_calc; [RTs, repmat(i, length(RTs), 1)]];
end

[pval_SFpq_r, r_SFpq_r]=permutation_pair_test_fast(SRTr_calc(:,1),SRTr_calc(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
nTemp = nTemp + length(SRTr_calc(:,1));

pqsingfracNR = SRTnr_calc;
pqsingfracR = SRTr_calc;

% Assemble data
searchNR(:,3)    = ones(1,length(searchNR))*1;
searchR(:,3)     = ones(1,length(searchR))*2;

singfracNR(:,3)  = ones(1,length(singfracNR))*3;
singfracR(:,3)   = ones(1,length(singfracR))*4;

pqsingfracNR(:,3)= ones(1,length(pqsingfracNR))*3;
pqsingfracR(:,3) = ones(1,length(pqsingfracR))*4;

% bootstrap errorbars
data_errBar = [searchNR];
corrDiff = @(z) corrMod(z(:,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_errBar},'alpha',0.001,'type','per');
errBar_searchNR3 = std(BOOTSTAT);

data_errBar = [searchR];
corrDiff = @(z) corrMod(z(:,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_errBar},'alpha',0.001,'type','per');
errBar_searchR3 = std(BOOTSTAT);

data_errBar = [singfracNR];
corrDiff = @(z) corrMod(z(:,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_errBar},'alpha',0.001,'type','per');
errBar_singfracNR3 = std(BOOTSTAT);

data_errBar = [singfracR];
corrDiff = @(z) corrMod(z(:,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_errBar},'alpha',0.001,'type','per');
errBar_singfracR3 = std(BOOTSTAT);

data_errBar = [pqsingfracNR];
corrDiff = @(z) corrMod(z(:,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_errBar},'alpha',0.001,'type','per');
errBar_pqsingfracNR3 = std(BOOTSTAT);

data_errBar = [pqsingfracR];
corrDiff = @(z) corrMod(z(:,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_errBar},'alpha',0.001,'type','per');
errBar_pqsingfracR3 = std(BOOTSTAT);

% bootErrBar3 = [errBar_searchNR3, errBar_searchR3, errBar_singfracNR3, errBar_singfracR3, errBar_pqsingfracNR3, errBar_pqsingfracR3];
bootErrBar3 = [errBar_pqsingfracNR3, errBar_pqsingfracR3, errBar_singfracNR3, errBar_singfracR3, errBar_searchNR3, errBar_searchR3];

%---------------------------------------
% Bootstrap of difference between SF and Search
% Note: corrMod is a wrapper to only return the corrVal from
% corr, defined at the end of the file

% Assemble data
searchNR(:,3)    = ones(1,length(searchNR))*1;
searchR(:,3)     = ones(1,length(searchR))*2;

singfracNR(:,3)  = ones(1,length(singfracNR))*3;
singfracR(:,3)   = ones(1,length(singfracR))*4;

pqsingfracNR(:,3)= ones(1,length(pqsingfracNR))*3;
pqsingfracR(:,3) = ones(1,length(pqsingfracR))*4;

% bootstrap sf vs search
data_sfVSsearchNR = [searchNR;singfracNR];
corrDiff = @(z) corrMod(z(z(:,3)==1,1:2)) - corrMod(z(z(:,3)==3,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_sfVSsearchNR},'alpha',0.01,'type','per');
data_sfVSsearchNR_CI = CI;

data_sfVSsearchR = [searchR;singfracR];
corrDiff = @(z) corrMod(z(z(:,3)==2,1:2)) - corrMod(z(z(:,3)==4,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_sfVSsearchR},'alpha',0.01,'type','per');
data_sfVSsearchR_CI = CI;

% bootstrap pqsf vs search

data_pqsfVSsearchNR = [searchNR;pqsingfracNR];
corrDiff = @(z) corrMod(z(z(:,3)==1,1:2)) - corrMod(z(z(:,3)==3,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_pqsfVSsearchNR},'alpha',0.01,'type','per');
data_pqsfVSsearchNR_CI = CI;

data_pqsfVSsearchR = [searchR;pqsingfracR];
corrDiff = @(z) corrMod(z(z(:,3)==2,1:2)) - corrMod(z(z(:,3)==4,1:2));
[CI,BOOTSTAT] = bootci(20000,{corrDiff, data_pqsfVSsearchR},'alpha',0.01,'type','per');
data_pqsfVSsearchR_CI = CI;

% Plot results
rOffset = -0.1;
nrOffset = -rOffset;
xs = [1+nrOffset, 2+rOffset, 3+nrOffset, 4+rOffset, 5+nrOffset, 6+rOffset];
astFontsize = 24;
textFontsize = 12;

mydata = [r_SFpq_nr, r_SFpq_r, r_SF_nr, r_SF_r, r_S_nr, r_S_r];
errBarLen = bootErrBar3;
colors = [nrewColor; rewColor; nrewColor; rewColor; nrewColor; rewColor];
for i = 1:6
    b = bar(xs(i), mydata(i)); hold on;
    set(b,'FaceColor',colors(i,:));
    errorbar(xs(i), mydata(i), errBarLen(i),  '.', 'Color', 'k');
end

set(gca, 'XTick', [1.5, 3.5, 5.5], 'XTickLabel', {'Environments Precued', 'Environment-Specific targets', 'Search'});

rref = [r_SFpq_nr, r_SFpq_r, r_SF_nr, r_SF_r, r_S_nr, r_S_r];
pref = [pval_SFpq_nr, pval_SFpq_r, pval_SF_nr, pval_SF_r, pval_S_nr, pval_S_r];
for i=1:6
    if pref(i) < 0.0001
        text(xs(i),rref(i)+0.01,'***','FontSize',astFontsize, 'HorizontalAlignment', 'center');
    elseif pref(i) < 0.05
        text(xs(i),rref(i)+0.01,'*','FontSize',astFontsize, 'HorizontalAlignment', 'center');
    else
        text(xs(i),rref(i)+0.01,'n.s.','FontSize',textFontsize, 'HorizontalAlignment', 'center');
    end
end

if all(data_sfVSsearchNR_CI>0) || all(data_sfVSsearchNR_CI<0) && all(data_sfVSsearchR_CI>0) || all(data_sfVSsearchR_CI<0)
    text(mean(xs(4:5)),.25,'***','FontSize',astFontsize, 'HorizontalAlignment', 'center');
else
    text(mean(xs(4:5)),.25,'n.s.','FontSize',textFontsize, 'HorizontalAlignment', 'center');
end
line([mean(xs(3:4)) mean(xs(5:6))], [.24 .24], 'Color', 'k');

if all(data_pqsfVSsearchNR_CI>0) || all(data_pqsfVSsearchNR_CI<0) && all(data_pqsfVSsearchR_CI>0) || all(data_pqsfVSsearchR_CI<0)
    text(mean(xs(2:5)),.27,'***','FontSize',astFontsize, 'HorizontalAlignment', 'center');
else
    text(mean(xs(2:5)),.27,'n.s.','FontSize',textFontsize, 'HorizontalAlignment', 'center');
end
line([mean(xs(1:2)) mean(xs(5:6))], [.26 .26], 'Color', 'k');

ylim([-0.0500    0.2900])
curr = ylim;
ylim([curr(1), curr(2)+0.02])
ylabel('Correlation (r)')

function retVal = rocarea3_area(a, b) % return only the first value from rocarea3 for area difference calculations
    temp = rocarea3(a,b);
    retVal = temp(1);
end

function retVal = rocWrapper(z)
    retVal = rocarea3(z(z(:,2)==1,1),z(z(:,2)==2,1));
end

function retVal = corrMod(a) 
    temp = corr(a, 'type', 'spearman');
    retVal = temp(1,2);
end