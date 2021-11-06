clear all; %close all hidden; clc;
addpath('functions'); 
load('Upload_storage_Fig1.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrewColor = [0.6235    0.6196    0.9961];
rewColor  = [0.9686    0.4627    0.5137];
meanColor = [1.0000    0.9961    0.3294];
rng(1);
%%
figure;
badRed = nrewColor;
goodBlue = rewColor;

nuDate = storage.Date(:,3)*1e4 + storage.Date(:,2)*1e2 + storage.Date(:,1)*1;
nuTime = storage.Date(:,4)*1e2 + storage.Date(:,5)*1;
nuID = nuDate*10 + storage.monkID';
uniqID = unique(nuID);
iter=0;
for session = (length(uniqID)-19):length(uniqID)
    iter = iter+1;
    % NO Rew
    tempIndex = (nuID'==uniqID(session))' & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'==0;
    st{1,iter} =  storage.searchDuration(tempIndex,1);
    
    % Rew
    tempIndex = (nuID'==uniqID(session))' & ((storage.finType'==0.2)|(storage.finType'==0.0)) & storage.rewDur'>0;
    st{2,iter} =  storage.searchDuration(tempIndex,1);
end

mdst = cellfun(@median, st);
mst = cellfun(@mean,st);
stdst = cellfun(@std,st);
sest = stdst./sqrt(cellfun(@length,st));

x = [1:20, 20:-1:1];
sespace = [mst+sest, fliplr(mst-sest)];

i = 1;
fill(x, sespace(i,:), badRed); hold on;
plot(1:size(mst(i,:),2), mst); hold on;

i = 2;
fill(x, sespace(i,:), goodBlue); hold on;
plot(1:size(mst(i,:),2), mst); hold on;
title('Search task')
ylabel('Search duration (s)');
xlabel('Last 20 Sessions');

%%
figure
iter=0;
for session = (length(uniqID)-19):length(uniqID)
    iter = iter+1;
    % NO Rew
    tempIndex = (nuID'==uniqID(session))' & ((storage.finType'==0.3)) & storage.rewDur'==0;
    stSF{1,iter} =  storage.searchDuration(tempIndex,1);
    
    % Rew
    tempIndex = (nuID'==uniqID(session))' & ((storage.finType'==0.3)) & storage.rewDur'>0;
    stSF{2,iter} =  storage.searchDuration(tempIndex,1);
end



mdstSF = cellfun(@median, stSF);
mstSF = cellfun(@mean,stSF);
stdstSF = cellfun(@std,stSF);
sestSF = stdstSF./sqrt(cellfun(@length,stSF));

x = [1:20, 20:-1:1];
sespaceSF = [mstSF+sestSF, fliplr(mstSF-sestSF)];

i = 1;
fill(x, sespaceSF(i,:), badRed); hold on;
plot(1:size(mstSF(i,:),2), mstSF); hold on;

i = 2;
fill(x, sespaceSF(i,:), goodBlue); hold on;
plot(1:size(mstSF(i,:),2), mstSF); hold on;
title('Single fractal task')
ylabel('Search duration (s)');
xlabel('Last 20 Sessions');
