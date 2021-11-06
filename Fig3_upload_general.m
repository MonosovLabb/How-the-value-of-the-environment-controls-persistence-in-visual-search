clear all;
addpath('functions'); 
addpath('MatlabToolbox-2.8');
load('Upload_storage_General.mat');

rng(1);

excelOut = 'tableTemp.xlsx';

nrewColor = [0.6235    0.6196    0.9961];
rewColor  = [0.9686    0.4627    0.5137];
sizeMeanMark = 20;
astFontsize = 24;
textFontsize = 12;

NumberofPerm=10000; %shuffle 10,000 times to get P-values for correlations


%% Panel A
fractalPosns = [22.5, 6;    112.5, 6;   202.5, 6;   292.5, 6    ;...    % First ring
    67.5, 15;   157.5, 15;  247.5, 15;  337.5, 15   ;...    % Second Ring
    22.5, 17;   112.5, 17;  202.5, 17;  292.5, 17   ];      % Third Ring
fractalPosns(fractalPosns(:,1)>180,1) = fractalPosns(fractalPosns(:,1)>180,1) -360;

for i = 1:5
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
end

% Initial Saccade Analysis
final = storage.targDir' - storage.firstSaccDir';
for i = 1:4
    % NO Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'==0 & ~storage.novelTrial' & storage.goodtrials');
    finalsorted{i,1,1} =  final(tempIndex);
    finalsorted{i,1,1} = mod(finalsorted{i,1,1}(~isnan(finalsorted{i,1,1})),360);
    finalsorted{i,1,1}(finalsorted{i,1,1}>180) = finalsorted{i,1,1}(finalsorted{i,1,1}>180)-360;
    
    % Rew
    tempIndex = logical(env(:,i+1) & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'>0 & ~storage.novelTrial' & storage.goodtrials');
    finalsorted{i+1,2,1} =  final(tempIndex);
    finalsorted{i+1,2,1} = mod(finalsorted{i+1,2,1}(~isnan(finalsorted{i+1,2,1})),360);
    finalsorted{i+1,2,1}(finalsorted{i+1,2,1}>180) = finalsorted{i+1,2,1}(finalsorted{i+1,2,1}>180)-360;
end

% Get plotting settings
maxAx = 0;
tempFig = figure('visible','off');
binCount = -175:10:175;
for i = 1:4
    histogram(finalsorted{i,1}, binCount, 'normalization', 'probability'); hold on;
    ax(i,1) = gca;
    maxAx = max(maxAx, ax(i,1).YLim(2));
    histogram(finalsorted{i+1,2}, binCount, 'normalization', 'probability'); hold on;
    ax(i+1,2) = gca;
    maxAx = max(maxAx, ax(i+1,2).YLim(2));
end
close(tempFig);

figure;
angOffset = 90;
labels = {'0%','25%','50%','75%','100%'};

rowCount = 2;
colCount = 20;
for i = 1:4
    subplot(rowCount, colCount, ((i*4)-3 + colCount):((i*4) + colCount));
    hold on;
    for iter = 1:4
        rectangle('Pos',[fractalPosns(iter,1)-5, 0, 10, maxAx/6],'FaceColor',[0.9686    0.9176    0.8431],'EdgeColor','none');
    end
    for iter = 5:8
        rectangle('Pos',[fractalPosns(iter,1)-5, 0, 10, maxAx/6],'FaceColor',[0.8784    0.7725    0.6078],'EdgeColor','none');
    end
    rectangle('Pos',[0-5, 0, 10, maxAx],'FaceColor',[0.9020    0.9020    0.9020],'EdgeColor','none');
    histogram(finalsorted{i,1}, binCount, 'normalization', 'probability','FaceColor',nrewColor,'EdgeColor','k'); hold on;
    ax(i,1) = gca;
    tempCurr = ax(i,1).YLim(2);
    ax(i,1).YLim = [0 maxAx];
    xticks(-90:90:90);
    xlim([-180 180]);
    yticks([0 0.8]);
    title(labels{i});
    if i==1
        ylabel('Proportion of initial saccades');
    end

    subplot(rowCount, colCount, ((i*4)+1):(((i+1)*4)));
    hold on;
    for iter = 1:4
        rectangle('Pos',[fractalPosns(iter,1)-5, 0, 10, maxAx/6],'FaceColor',[0.9686    0.9176    0.8431],'EdgeColor','none');
    end
    for iter = 5:8
        rectangle('Pos',[fractalPosns(iter,1)-5, 0, 10, maxAx/6],'FaceColor',[0.8784    0.7725    0.6078],'EdgeColor','none');
    end
    rectangle('Pos',[0-5, 0, 10, maxAx],'FaceColor',[0.9020    0.9020    0.9020],'EdgeColor','none');
    histogram(finalsorted{i+1,2}, binCount, 'normalization', 'probability','FaceColor',rewColor,'EdgeColor','k'); hold on;
    ax(i+1,2) = gca;
    tempCurr = ax(i+1,2).YLim(2);
    ax(i+1,2).YLim = [0 maxAx];
    xticks(-90:90:90);
    xlim([-180 180]);
    yticks([0 0.8]);
    title(labels{i+1});
    if i==1
        ylabel('Proportion of initial saccades');
    end
end
sgtitle({'The majority of initial saccades are directed toward the target fractal,' 'especially for the Reward target in low-reward environments'})
drawnow(); scale(gcf, 'hscale',3)
%% Panel B
% Proportion direct 1.5 deg.
figure;

propFirstSucc = {};

for i = 1:5
    env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i))';% & nuID==uniqID(session);
    % NO Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'==0 & ~storage.novelTrial' & storage.goodtrials');
    temp =  storage.firstSacc2Targ3(1,tempIndex);
    temp =  temp(~isnan(temp));
    propFirstSucc{i,1} = temp;
    
    % Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'>0 & ~storage.novelTrial' & storage.goodtrials');
    temp = storage.firstSacc2Targ3(1,tempIndex);
    temp = temp(~isnan(temp));
    propFirstSucc{i,2} = temp;
end

mst = cellfun(@mean, propFirstSucc);
for iter = 1:4
    sest(iter,1) = std(bootstrp(20000, @mean, propFirstSucc{iter,1}));
    sest(iter+1,2) = std(bootstrp(20000, @mean, propFirstSucc{iter+1,2}));
end
propSucc_NR = [];
propSucc_R =  [];

for iter = 1:4
    propSucc_NR = [propSucc_NR; [propFirstSucc{iter,1}', iter*ones(size(propFirstSucc{iter,1}'))]];
    propSucc_R = [propSucc_R; [propFirstSucc{iter+1,2}', (iter+1)*ones(size(propFirstSucc{iter+1,2}'))]];
end

[pval_r, r_r]=permutation_pair_test_fast(propSucc_R(:,1),propSucc_R(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
[pval_nr, r_nr]=permutation_pair_test_fast(propSucc_NR(:,1),propSucc_NR(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear

% Mean Search time & STD
rOffset = -0.1;
nrOffset = -rOffset;
xs = [1+nrOffset, 2+rOffset, 3+nrOffset, 4+rOffset, 5+nrOffset, 6+rOffset, 7+nrOffset, 8+rOffset, 9+nrOffset, 10+rOffset];

hold on;
for i = 1:4
    % NoRew
    a(i,1) = bar(xs((2*i)-1), mst(i,1,1));
    a(i,1).FaceColor = nrewColor;
    errorbar(mst(i,1,1), sest(i,1,1), '.', 'Color','k', 'XData', xs((2*i)-1), 'LineWidth', 3);
    % Rew
    a(i,2) = bar(xs(2*i)+2, mst(i+1,2,1));
    a(i,2).FaceColor = rewColor;
    errorbar(mst(i+1,2,1), sest(i+1,2,1), '.', 'Color','k', 'XData', xs(2*i)+2, 'LineWidth', 3);
end

Nrew = errorbar(mst(1:4,1,1), [0 0 0 0], '-', 'Color',nrewColor, 'XData', xs((2*(1:4))-1), 'LineWidth', 3);
Rrew = errorbar(mst((1:4)+1,2,1), [0 0 0 0], '-', 'Color',rewColor, 'XData', xs(2*(1:4))+2, 'LineWidth', 3);

text(1,0.9,['p = ' mat2str(pval_nr,3)], 'Color', nrewColor);
text(1,0.8,['r = ' mat2str(r_nr,3)], 'Color', nrewColor);
text(5,0.9,['p = ' mat2str(pval_r,3)], 'Color', rewColor);
text(5,0.8,['r = ' mat2str(r_r,3)], 'Color', rewColor);
xlim([0.5 10.5])
currY = ylim;
ylim([0 1])
xlabel('Environment reward probability')
ylabel('Proportion of initial saccades direct to the target fractal')
title({'Direct saccades to the target' 'increase with its reward value and' 'decrease with reward probability'})
set(gca, 'XTick', [1.5, 3.5, 5.5, 7.5, 9.5, 11.5], 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
hold off

%% Panel C
%% First saccade 2 target persistent search time - Full Time with timeout trials
figure;

propFirstSucc = {};
SRT_calc_NR = [];
SRT_calc_R  = [];

rewFill = rewColor;
nRewFill = nrewColor;
linew = 1;

xOffset = 0.1;
noRewXs = 2*(1:4)-1+xOffset;
rewXs = 2*((1:4))+2 - xOffset;

for i = 1:5
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    % NO Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'==0 & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
    propFirstSucc{i,1,1} =  storage.firstSacc2TargPersSearchDur(tempIndex);
    propFirstSucc{i,1,1} =  propFirstSucc{i,1,1}(~isnan(propFirstSucc{i,1,1}));
    SRT_calc_NR = [SRT_calc_NR; [propFirstSucc{i,1,1}', repmat(i, size(propFirstSucc{i,1,1}',1), 1)]];
    % Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'>0 & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
    propFirstSucc{i,2,1} =  storage.firstSacc2TargPersSearchDur(tempIndex);
    propFirstSucc{i,2,1} =  propFirstSucc{i,2,1}(~isnan(propFirstSucc{i,2,1}));
    SRT_calc_R =  [SRT_calc_R; [propFirstSucc{i,2,1}', repmat(i, size(propFirstSucc{i,2,1}',1), 1)]];
end

% propFirstSucc
mpropFirstSucc = cellfun(@mean,propFirstSucc);
stdpropFirstSucc = cellfun(@std,propFirstSucc);
sepropFirstSucc = stdpropFirstSucc./sqrt(cellfun(@length,propFirstSucc));

% Significance analyses
for i = 1:3
    pstor(i,1) = ranksum(propFirstSucc{i,1,1}, propFirstSucc{i+1,1,1});
    pstor(i+1,2) = ranksum(propFirstSucc{i+1,2,1}, propFirstSucc{i+2,2,1});
end
[pvalNR, rNR] = permutation_pair_test_fast(SRT_calc_NR(:,1),SRT_calc_NR(:,2),NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear
[pvalR,  rR ] = permutation_pair_test_fast(SRT_calc_R(:,1),SRT_calc_R(:,2),NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear

% Mean Search time & propFirstSuccD
hold on;
for i = 1:4
    % NoRew
    iosr.statistics.boxPlot(propFirstSucc{i,1}', 'showviolin', false, 'boxcolor', nRewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', noRewXs(i), 'boxWidth', 0.8, 'showMean', true, 'meanColor',  'b', 'meanMarker', '.', 'meanSize', sizeMeanMark); hold on;
    
    % Rew
    iosr.statistics.boxPlot(propFirstSucc{i+1,2}', 'showviolin', false, 'boxcolor', rewFill , 'notch', false, 'linewidth', linew, 'showOutliers', false, 'medianColor', 'k', 'lineColor', 'k', 'x', rewXs(i), 'boxWidth', 0.8, 'showMean', true, 'meanColor',  'r', 'meanMarker', '.', 'meanSize', sizeMeanMark ); hold on;
end

xlim([0.5 10.5])

xlabel('Environment reward probability')
ylabel({'Search duration after' 'gazing at target (s)'});
yticks(0:5);
yticklabels({'0' , '',  '',  '',  '', '5'});
title({'Persistent search after direct saccades' 'to the No Reward target' 'increases with reward probability'});
set(gca, 'XTick', [1.5, 3.5, 5.5, 7.5, 9.5, 11.5], 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});

clear net sing
[net, sing] = chartCalcs(propFirstSucc, NumberofPerm);
textPrintMeans(propFirstSucc, net, sing, noRewXs, rewXs, nan, nrewColor, rewColor);

hold off

%% Panel D
clearvars -except rewColor nrewColor NumberofPerm
rng(1);
load('Upload_storage_Fig3Novel.mat')

figure

noRewXs = 2*(1:4)-0.75;
rewXs = 2*((1:4))+1.75;
allXs = [1.5, 3.5, 5.5, 7.5, 9.5, 11.5];
rOffset = -0.1;
nrOffset = -rOffset;
xs = [1+nrOffset, 2+rOffset, 3+nrOffset, 4+rOffset, 5+nrOffset, 6+rOffset, 7+nrOffset, 8+rOffset, 9+nrOffset, 10+rOffset];

objFixed = {};
for i = 1:5
    env(:,i) = storage.envID==(0+i) | storage.envID==(6+i);
    % NO Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)|(storage.finType'==3.0)) & storage.rewDur'==0 & storage.novelTrial' & storage.goodtrials');
    objFixed{i,1,1} =  storage.durFixNovel(tempIndex)./(storage.searchDuration(tempIndex,1))';
    objFixed{i,1,1} =  objFixed{i,1,1}(~isnan(objFixed{i,1,1}));
    % Rew
    tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0)|(storage.finType'==3.0)) & storage.rewDur'>0 & storage.novelTrial' & storage.goodtrials');
    objFixed{i,2,1} =  storage.durFixNovel(tempIndex)./(storage.searchDuration(tempIndex,1))';
    objFixed{i,2,1} =  objFixed{i,2,1}(~isnan(objFixed{i,2,1}));
end

mobjFixed = cellfun(@mean,objFixed);
stdobjFixed = cellfun(@std,objFixed);
seobjFixed = stdobjFixed./sqrt(cellfun(@length,objFixed));

hold on;
Nrew = errorbar(mobjFixed(1:4,1,1), [0 0 0 0], '-', 'Color',nrewColor, 'XData', xs((2*(1:4))-1), 'LineWidth', 3);
for i = 1:4
    % NoRew
    a(i,1) = bar(xs((2*i)-1), mobjFixed(i,1,1));
    a(i,1).FaceColor = nrewColor;
    errorbar(mobjFixed(i,1,1), seobjFixed(i,1,1), '.', 'Color','k', 'XData', xs((2*i)-1), 'LineWidth', 1.5);
end
Rrew = errorbar(mobjFixed((1:4)+1,2,1), [0 0 0 0], '-', 'Color',rewColor, 'XData', xs(2*(1:4))+2, 'LineWidth', 3);
for i = 1:4
    % Rew
    a(i,2) = bar(xs(2*i)+2, mobjFixed(i+1,2,1));
    a(i,2).FaceColor = rewColor;
    errorbar(mobjFixed(i+1,2,1), seobjFixed(i+1,2,1), '.', 'Color','k', 'XData', xs(2*i)+2, 'LineWidth', 1.5);
end

xlim([0.5 10.5])
currY = ylim;
xlabel('Environment reward probability')

yticks([0 0.25]);
ylabel({'Proportion of search duration' 'with gaze on novel fractal'});

title({'Gaze at a distracting novel' 'fractal on No Reward target trials' 'decreases with reward probability'});
set(gca, 'XTick', [1.5, 3.5, 5.5, 7.5, 9.5, 11.5], 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});

clear net sing
[net, sing] = chartCalcs(objFixed, NumberofPerm);
textPrintMeans(objFixed, net, sing, noRewXs, rewXs, 0.1, nrewColor, rewColor);

hold off

