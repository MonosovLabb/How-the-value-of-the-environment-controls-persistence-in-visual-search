clear all;
addpath('functions'); 
addpath('MatlabToolbox-2.8');
rng(1);
load('Upload_storage_General.mat');

% Switch IDs to simplify plotting sequence
storage.monkID(storage.monkID==1) = 3;
storage.monkID(storage.monkID==2) = 1;
storage.monkID(storage.monkID==3) = 2;

figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrewColor = [0.6235    0.6196    0.9961];
rewColor  = [0.9686    0.4627    0.5137];
meanColor = [1.0000    0.9961    0.3294];
rowCount = 7;
colCount = 20;
widScaler = 4;
sizeMeanMark = 10;

NumberofPerm=10000; %need to shuffle 10,000 times to get real P-values for correlations


%%%%plotparameters
BinWidth=0.025; % histogram bin
XLim=2;
JitterVar=0.3; % jitter for scatter plot
SizeofScat=1;

animNames = {'Animal B' 'Animal M'};

%% Target Reaction Time
for monkParam = 1:2
    Monkid = monkParam;
    noRewXs = 2*(1:4)-0.75;
    rewXs = 2*((1:4))+1.75;
    allXs = [1.5, 3.5, 5.5, 7.5, 9.5, 11.5];
    propFirstSucc = {};
    SRT_calc_NR = [];
    SRT_calc_R  = [];
    SRT_calc_Both = [];
    env = [];

    linew = 1;
    
    xOffset = 0.1;
    noRewXs = 2*(1:4)-1+xOffset;
    rewXs = 2*((1:4))+2 - xOffset;
    
    for i = 1:5
        switch monkParam
            case 1
                env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & storage.monkID==1;
            case 2
                env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & storage.monkID==2;
        end
        
        % NO Rew
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'==0 & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
        tempIndex2 = tempIndex & ~isnan(storage.targRT)';
        propFirstSucc{i,1,1} =  storage.targRT(tempIndex2);
        propFirstSucc{i,1,1} =  propFirstSucc{i,1,1}(~isnan(propFirstSucc{i,1,1}));
        SRT_calc_NR = [SRT_calc_NR; [propFirstSucc{i,1,1}', repmat(i, size(propFirstSucc{i,1,1}',1),1)]];
        
        % Rew
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'>0 & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
        tempIndex2 = tempIndex & ~isnan(storage.targRT)';
        propFirstSucc{i,2,1} =  storage.targRT(tempIndex2);
        propFirstSucc{i,2,1} =  propFirstSucc{i,2,1}(~isnan(propFirstSucc{i,2,1}));
        SRT_calc_R =  [SRT_calc_R; [propFirstSucc{i,2,1}', repmat(i, size(propFirstSucc{i,2,1}',1), 1)]];
        
        % Both
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.firstSacc2Targ3(1,:)'==1 & ~storage.novelTrial' & storage.goodtrials');
        tempIndex2 = tempIndex & ~isnan(storage.targRT)';
        propFirstSucc{i,3,1} =  storage.targRT(tempIndex2);
        propFirstSucc{i,3,1} =  propFirstSucc{i,3,1}(~isnan(propFirstSucc{i,3,1}));
        SRT_calc_Both =  [SRT_calc_Both; [propFirstSucc{i,3,1}', repmat(i, size(propFirstSucc{i,3,1}',1), 1)]];
    end
    
    mst = cellfun(@mean,propFirstSucc);
    stdst = cellfun(@std,propFirstSucc);
    sest = stdst./sqrt(cellfun(@length,propFirstSucc));
    
    % Significance analyses
    pstor = [];
    for i = 1:3
        pstor(i,1) = ranksum(propFirstSucc{i,1,1}, propFirstSucc{i+1,1,1});
        pstor(i+1,2) = ranksum(propFirstSucc{i+1,2,1}, propFirstSucc{i+2,2,1});
    end
    for i = 1:4
        pstor(i,3) = ranksum(propFirstSucc{i,3,1}, propFirstSucc{i+1,3,1});
    end
    
    pval = ranksum(propFirstSucc{4,3}, propFirstSucc{5,3});
    disp(['For '  animNames{monkParam} ' the ranksum diff between 75% and 100% reaction time p=' num2str(pval)]);
    
    [pval_nr, r_nr] = permutation_pair_test_fast(SRT_calc_NR(:,1),SRT_calc_NR(:,2),NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear
    [pval_r,  r_r ] = permutation_pair_test_fast(SRT_calc_R(:,1),SRT_calc_R(:,2),NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear
    [pval_both,  r_both] = permutation_pair_test_fast(SRT_calc_Both(:,1),SRT_calc_Both(:,2),NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear
    
    clear net sing
    [net, sing] = chartCalcs(propFirstSucc, NumberofPerm);
    
    subplot(4, 2, monkParam+6);
    hold on;
    % NoRew
    errorbar(mst(1:4,1,1), sest(1:4,1,1), '-', 'Color',nrewColor, 'XData', 1:4, 'LineWidth', 3);
    % Rew
    errorbar(mst(2:5,2,1), sest(2:5,2,1), '-', 'Color',rewColor, 'XData', 2:5, 'LineWidth', 3);
    base = ylim();
    text(1,0.9*range(ylim)+base(1),['p_N_R = ' mat2str(pval_nr,3)]);
    text(1,0.8*range(ylim)+base(1),['r_N_R = ' mat2str(r_nr,3)]);
    text(3,0.9*range(ylim)+base(1),['p_R = ' mat2str(pval_r,3)]);
    text(3,0.8*range(ylim)+base(1),['r_R = ' mat2str(r_r,3)]);
    xlim([0.5 5.5])
    currY = ylim;
    xlabel('Environment Reward Probability')
    ylabel('Direct Saccade Response Time (s)')
    switch Monkid
        case 2
            name='Animal M'
        case 1
            name='Animal B'
    end
    title([name, ': Target Saccade Reaction Time']);
    set(gca, 'XTick', 1:5, 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
    if Monkid==1
        ylim([0.185 0.225])
    else
        ylim([0.22 0.34])
    end
    hold off;
    
    subplot(4, 2, monkParam+2);
    hold on;
    % Both
    errorbar(mst(1:5,3,1), sest(1:5,3,1), '-', 'Color','k', 'XData', 1:5, 'LineWidth', 3);
    base = ylim();
    text(5,0.9*range(ylim)+base(1),['p_B_o_t_h = ' mat2str(pval_both,3)]);
    text(5,0.8*range(ylim)+base(1),['r_B_o_t_h = ' mat2str(r_both,3)]);
    hold off;
    
    xlim([0.5 5.5])
    currY = ylim;
    xlabel('Environment Reward Probability')
    ylabel('Direct Saccade Response Time (s)')
    switch Monkid
        case 2
            name='Animal M'
        case 1
            name='Animal B'
    end
    title([name, ': Target Saccade Reaction Time']);
    set(gca, 'XTick', 1:5, 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
    hold off
    if Monkid==1
        ylim([0.185 0.225])
    else
        ylim([0.22 0.34])
    end
end

%% Proportion direct 1.5 deg.
for monkParam = 1:2
    Monkid = monkParam;
    
    badRed = [0.9451 0.2627 0.2471];
    goodBlue = [0.4392 0.7176 0.7294];
    noRewXs = 2*(1:4)-0.75;
    rewXs = 2*((1:4))+1.75;
    allXs = [1.5, 3.5, 5.5, 7.5, 9.5, 11.5];
    
    propFirstSucc = {};
    SRT_calc_NR = [];
    SRT_calc_R  = [];
    env = [];
    
    for i = 1:5
        switch monkParam
            case 1
                env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & storage.monkID==1;
            case 2
                env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i)) & storage.monkID==2;
        end
        %     env(:,i) = (storage.envID==(0+i) | storage.envID==(6+i))';% & nuID==uniqID(session);
        % NO Rew
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'==0 & ~storage.novelTrial' & storage.goodtrials');
        tempIndex2 = tempIndex & ~isnan(storage.firstSacc2Targ3(1,:)');
        temp =  storage.firstSacc2Targ3(1,tempIndex2);
        propFirstSucc{i,1} = temp;
        
        % Rew
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & storage.rewDur'>0 & ~storage.novelTrial' & storage.goodtrials');
        tempIndex2 = tempIndex & ~isnan(storage.firstSacc2Targ3(1,:)');
        temp = storage.firstSacc2Targ3(1,tempIndex2);
        propFirstSucc{i,2} = temp;
        
        % Both
        tempIndex = logical(env(:,i) & ((storage.finType'==0.2)|(storage.finType'==0.0|storage.finType'==3.0)) & ~storage.novelTrial' & storage.goodtrials');
        tempIndex2 = tempIndex & ~isnan(storage.firstSacc2Targ3(1,:)');
        temp = storage.firstSacc2Targ3(1,tempIndex2);
        propFirstSucc{i,3} = temp;
    end

    % Errorbars
    mst = cellfun(@mean, propFirstSucc);
    sest = [];
    for iter = 1:4
        sest(iter,1) = std(bootstrp(20000, @mean, propFirstSucc{iter,1}));
        sest(iter+1,2) = std(bootstrp(20000, @mean, propFirstSucc{iter+1,2}));
    end
    for iter = 1:5
        sest(iter,3) = std(bootstrp(20000, @mean, [propFirstSucc{iter,1}, propFirstSucc{iter,2}]));
    end
    
    % Significance analyses
    propSucc_NR = [];
    propSucc_R =  [];
    propSucc_Both =  [];
    for iter = 1:4
        propSucc_NR = [propSucc_NR; [propFirstSucc{iter,1}', iter*ones(size(propFirstSucc{iter,1}'))]];
        propSucc_R = [propSucc_R; [propFirstSucc{iter+1,2}', (iter+1)*ones(size(propFirstSucc{iter+1,2}'))]];
    end
    for iter = 1:5
        propSucc_Both = [propSucc_Both; [propFirstSucc{iter,3}', (iter+1)*ones(size(propFirstSucc{iter,3}'))]];
    end
    pval = ranksum(propFirstSucc{4,3}, propFirstSucc{5,3});
    disp(['For '  animNames{monkParam} ' the ranksum diff between 75% and 100% proportion direct p=' num2str(pval)]);
    [pval_r, r_r]=permutation_pair_test_fast(propSucc_R(:,1),propSucc_R(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
    [pval_nr, r_nr]=permutation_pair_test_fast(propSucc_NR(:,1),propSucc_NR(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
    [pval_both, r_both]=permutation_pair_test_fast(propSucc_Both(:,1),propSucc_Both(:,2),NumberofPerm,'rankcorr'); %corr for non parametric rank based; corr is linear
    
    pstor = [];
    for i = 1:3
        pstor(i,1) = ranksum(propFirstSucc{i,1,1}, propFirstSucc{i+1,1,1});
        pstor(i+1,2) = ranksum(propFirstSucc{i+1,2,1}, propFirstSucc{i+2,2,1});
    end
    for i = 1:4
        pstor(i,3) = ranksum(propFirstSucc{i,3,1}, propFirstSucc{i+1,3,1});
    end
    
    % Mean Search time & STD
    rOffset = -0.1;
    nrOffset = -rOffset;
    xs = [1+nrOffset, 2+rOffset, 3+nrOffset, 4+rOffset, 5+nrOffset, 6+rOffset, 7+nrOffset, 8+rOffset, 9+nrOffset, 10+rOffset];
    
    subplot(4, 2, monkParam+4);
    hold on;
    % NoRew
    errorbar(mst(1:4,1,1), sest(1:4,1,1), '-', 'Color',nrewColor, 'XData', 1:4, 'LineWidth', 3);
    % Rew
    errorbar(mst(2:5,2,1), sest(2:5,2,1), '-', 'Color',rewColor, 'XData', 2:5, 'LineWidth', 3);
    text(1,0.9,['p_N_R = ' mat2str(pval_nr,3)]);
    text(1,0.8,['r_N_R = ' mat2str(r_nr,3)]);
    text(3,0.9,['p_R = ' mat2str(pval_r,3)]);
    text(3,0.8,['r_R = ' mat2str(r_r,3)]);
    xlim([0.5 5.5])
    currY = ylim;
    ylim([0 1])
    xlabel('Environment Reward Probability')
    ylabel('Direct Saccade Probability')
    switch Monkid
        case 2
            name='Animal M'
        case 1
            name='Animal B'
    end
    title([name, ': Proportion Saccades Directly to Target 1.5 deg.'])
    set(gca, 'XTick', 1:5, 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
    ylim([0.3 0.9])
    hold off;
    
    subplot(4, 2, monkParam+0);
    hold on;
    % Both
    errorbar(mst(1:5,3,1), sest(1:5,3,1), '-', 'Color','k', 'XData', 1:5, 'LineWidth', 3);
    text(5,0.9,['p_B_o_t_h = ' mat2str(pval_both,3)]);
    text(5,0.8,['r_B_o_t_h = ' mat2str(r_both,3)]);
    hold off;

    xlim([0.5 5.5])
    currY = ylim;
    ylim([0 1])
    xlabel('Environment Reward Probability')
    ylabel('Direct Saccade Probability')
    switch Monkid
        case 2
            name='Animal M'
        case 1
            name='Animal B'
    end
    title([name, ': Proportion Saccades Directly to Target 1.5 deg.'])
    set(gca, 'XTick', 1:5, 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
    ylim([0.3 0.9])
end