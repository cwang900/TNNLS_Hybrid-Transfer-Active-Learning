
clear
clc

%% Run code to generate results

GenerateResultsTrigRandom % This provides RandomTrig.mat
GenerateResultsTrigAL % This provides ALTrig.mat

% ALso please run the python code SimulationTrig for NN

%% Mean RMSE of all methods for random batches

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

% load results for NN from corresponding folder
RMSENN = readtable("...\RandomTrig_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = mean(RMSENN(RMSENN(:,1)~=0,Ind))';
load('RandomTrig.mat')
eOffline = mean(RMSEOffline(RMSEOffline(:,1)~=0,Ind))';
eNT = mean(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind))';
eST = mean(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind))';
eLMC = mean(RMSELMC(RMSELMC(:,1)~=0,Ind))';
eProp = mean(RMSE(RMSE(:,1)~=0,Ind))';

figure
t = tiledlayout(2,2);
nexttile % Mean RMSE of all methods for random batches
hold on
p1 = plot(eNT,...
    'LineStyle','-','Color',[255 0 255 255]/255,...
    'Marker','s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(eST,...
    'LineStyle','-','Color',[255 0 0 255]/255,...
    'Marker','d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(eLMC,...
    'LineStyle','-','Color',[197 90 17 255]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(eNN,...
    'LineStyle','-','Color',[153 67 96 255]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(eOffline,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0.2 1])
xticks([1 5:5:30])
yticks(0.2:0.2:1)
xticklabels('')
xlabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(a) RMSE for Random Batches','FontWeight','bold')
box on
grid on
legend([p1,p2,p3,p4,p5,p6],{'SPGP','S2S','LMC','NN','Offline','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)


%% Mean RMSE of all methods for Active batches

% load results for NN from corresponding folder
RMSENN = readtable("...\Active3D2S_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = mean(RMSENN(RMSENN(:,1)~=0,Ind))';
load('RandomTrig.mat')
eOffline = mean(RMSEOffline(RMSEOffline(:,1)~=0,Ind))';
load('ALTrig.mat')
eNT = mean(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind))';
eST = mean(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind))';
eLMC = mean(RMSELMC(RMSELMC(:,1)~=0,Ind))';
eProp = mean(RMSE(RMSE(:,1)~=0,Ind))';

nexttile % 
hold on
p1 = plot(eNT,...
    'LineStyle','-','Color',[255 0 255 255]/255,...
    'Marker','s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(eST,...
    'LineStyle','-','Color',[255 0 0 255]/255,...
    'Marker','d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(eLMC,...
    'LineStyle','-','Color',[197 90 17 255]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(eNN,...
    'LineStyle','-','Color',[153 67 96 255]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(eOffline,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0.2 1])
xticks([1 5:5:30])
yticks(0.2:0.2:1)
% xticklabels('')
xlabel('Round','FontSize',12,'FontWeight','bold')
ylabel('RMSE','FontSize',12,'FontWeight','bold')
title('(b) RMSE for Active Batches','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

%% Average run time of all methods for Active batches

tEndNN = readtable("...\ActiveTrig_RMSE.csv");
tEndNN = tEndNN{2:27,:};
tEndNN = mean(tEndNN(tEndNN(:,2)~=0,:),2);
load('RandomTrig.mat')
tEndOffline = mean(tEnd(tEnd(:,2)~=0,Ind));
load('ALTrig.mat')
tEndNonTransfer= mean(tEndNonTransfer(tEndNonTransfer(:,2)~=0,Ind));
tEndSGPTransfer = mean(tEndSGPTransfer(tEndSGPTransfer(:,2)~=0,Ind));
tEndLMC = mean(tEndLMC(tEndLMC(:,2)~=0,Ind));
tEndProp = mean(tEndProp(tEndProp(:,2)~=0,Ind));

nexttile % 
hold on
p1 = plot(tEndNonTransfer,...
    'm-s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(tEndSGPTransfer,...
    'r-d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(tEndLMC,...
    'LineStyle','-','Color',[197 90 17]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(tEndNN,...
    'LineStyle','-','Color',[153 67 96]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(tEndOffline,...
    'k-o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(tEndProp,...
    'b->','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 220])
xticks([1 5:5:30])
yticks([0:50:220])
xticklabels('')
% xlabel('Round','FontSize',12)
ylabel('Time (s)','FontSize',12)
box on
grid on
title('(c) Computational Time','FontSize',12,'FontWeight','bold')
% legend([p1,p2,p3,p4,p5,p6],{'SPGP','S2S','LMC','NN','Offline','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)%,'Offline'

%% Zoom in plot of computational time for all online methods

tEndNN = readtable("...\ActiveTrig_Time.csv");
tEndNN = tEndNN{2:27,:};
tEndNN = mean(tEndNN(tEndNN(:,2)~=0,:),2);
load('RandomTrig.mat')
tEndOffline = mean(tEnd(tEnd(:,2)~=0,Ind));
load('ALTrig.mat')
tEndNonTransfer= mean(tEndNonTransfer(tEndNonTransfer(:,2)~=0,Ind));
tEndSGPTransfer = mean(tEndSGPTransfer(tEndSGPTransfer(:,2)~=0,Ind));
tEndLMC = mean(tEndLMC(tEndLMC(:,2)~=0,Ind));
tEndProp = mean(tEndProp(tEndProp(:,2)~=0,Ind));

nexttile
hold on
p1 = plot(tEndNonTransfer,...
    'm-s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(tEndSGPTransfer,...
    'r-d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(tEndLMC,...
    'LineStyle','-','Color',[197 90 17]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(tEndNN,...
    'LineStyle','-','Color',[153 67 96]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(tEndProp,...
    'b->','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 8])
xticks([1 5:5:30])
yticks([0:2:8])
% xticklabels('')
xlabel('Round','FontSize',12)
ylabel('Time (s)','FontSize',12)
box on
grid on
title('','FontSize',12,'FontWeight','bold')
% legend([p1,p2,p3,p4,p5,p6],{'SPGP','S2S','LMC','NN','Offline','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)%,'Offline'

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% Paired t-test

% RMSE Random batches
RMSENN = readtable("...\RandomTrig_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = RMSENN(RMSENN(:,1)~=0,:);
eNN = eNN(:,end);
load('RandomTrig.mat')
eOffline = RMSEOffline(RMSEOffline(:,1)~=0,Ind);
eOffline = eOffline(:,end);
eNT = RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = RMSELMC(RMSELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = RMSE(RMSE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST eLMC eNN eOffline eProp];

pVRandomRMSE = zeros(1,5);
for i = 1:5
    [~,pVRandomRMSE(i)] = ttest(eAll(:,i),eAll(:,end));
end

% RMSE Active batches
RMSENN = readtable("...\ALTrig_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = RMSENN(RMSENN(:,1)~=0,:);
eNN = eNN(:,end);
load('RandomTrig.mat')
eOffline = RMSEOffline(RMSEOffline(:,1)~=0,Ind);
eOffline = eOffline(:,end);
load('ALTrig.mat')
eNT = RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = RMSELMC(RMSELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = RMSE(RMSE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST eLMC eNN eOffline eProp];

pVActiveRMSE = zeros(1,5);
for i = 1:5
    [~,pVActiveRMSE(i)] = ttest(eAll(:,i),eAll(:,end));
end


% MAPE Random batches
MAPENN = readtable("...\RandomTrig_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = MAPENN(MAPENN(:,1)~=0,:);
eNN = eNN(:,end);
load('RandomTrig.mat')
eOffline = MAPEOffline(MAPEOffline(:,1)~=0,Ind);
eOffline = eOffline(:,end);
eNT = MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = MAPELMC(MAPELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = MAPE(MAPE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST eLMC eNN eOffline eProp];

pVRandomMAPE = zeros(1,5);
for i = 1:5
    [~,pVRandomMAPE(i)] = ttest(eAll(:,i),eAll(:,end));
end

% MAPE Active batches
MAPENN = readtable("...\ActiveTrig_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = MAPENN(MAPENN(:,1)~=0,:);
eNN = eNN(:,end);
load('RandomTrig.mat')
eOffline = MAPEOffline(RMSEOffline(:,1)~=0,Ind);
eOffline = eOffline(:,end);
load('ALTrig.mat')
eNT = MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = MAPELMC(MAPELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = MAPE(MAPE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST eLMC eNN eOffline eProp];

pVActiveMAPE = zeros(1,5);
for i = 1:5
    [~,pVActiveMAPE(i)] = ttest(eAll(:,i),eAll(:,end));
end


%% Quantiles of RMSE of all methods for random batches

for k = 1

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

strRandomRMSE = cell(1,length(pVActiveRMSE));
for i = 1:length(pVActiveRMSE)
    if pVRandomRMSE(i)<0.001
        strRandomRMSE{i} = sprintf('%.3e', pVRandomRMSE(i));
    else
        strRandomRMSE{i} = sprintf('%.3f', pVRandomRMSE(i));
    end
end

RMSENN = readtable("...\RandomTrig_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = median(RMSENN(RMSENN(:,1)~=0,Ind))';
load('RandomTrig.mat')
eOffline = median(RMSEOffline(RMSEOffline(:,1)~=0,Ind))';
fOffline = [quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(4))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(5))'];
eNT = median(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind))';
eST = median(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind))';
eLMC = median(RMSELMC(RMSELMC(:,1)~=0,Ind))';
eProp = median(RMSE(RMSE(:,1)~=0,Ind))';

fNN = [quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(5))')];
fNT = [quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(5))')];
fST = [quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(5))')];
fLMC = [quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(5))')];
fProp = [quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(4))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(5))'];

figure
t = tiledlayout(4,2);
nexttile%([2 2]) % Mean RMSE of all methods for random batches
hold on
p1 = plot(eNT,...
    'LineStyle','-','Color',[255 0 255 255]/255,...
    'Marker','s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(eST,...
    'LineStyle','-','Color',[255 0 0 255]/255,...
    'Marker','d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(eLMC,...
    'LineStyle','-','Color',[197 90 17 255]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(eNN,...
    'LineStyle','-','Color',[153 67 96 255]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(eOffline,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 1.5])
xticks([1 5:5:30])
yticks(0:0.5:1.5)
xticklabels('')
xlabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(a) Median RMSE for Random Batches','FontWeight','bold')
box on
grid on
% legend([p1,p2,p3,p4,p5,p6],{'SPGP','S2S','LMC','NN','Offline','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for SPGP
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNT(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNT(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('','FontSize',12)
title(['(b) SPGP (p=' strRandomRMSE{1} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for S2S
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fST(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fST(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('RMSE','FontSize',12)
title(['(c) S2S (p=' strRandomRMSE{2} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for LMC
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{i},'Color',[197 90 17 255]/255,...
            'LineWidth',1.5);
    else
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{end},'Color',[197 90 17 255]/255 ...
            ,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('','FontSize',12)
title(['(d) LMC (p=' strRandomRMSE{3} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for NN
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNN(:,i),...
            'LineStyle',LineStyle{i},'Color',[153 67 96 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNN(:,i),...
            'LineStyle',LineStyle{end},'Color',[153 67 96 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 5])
xticks([1 5:5:30])
yticks(0:1:5)
xticklabels('')
xlabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title(['(e) NN (p=' strRandomRMSE{4} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for Offline
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('','FontSize',12)
ylabel('','FontSize',12)
title(['(f) Offline (p=' strRandomRMSE{5} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for Prop
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fProp(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fProp(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xlabel('Round','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(g) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

end

%% Quantiles of RMSE of all methods for Active batches

for k = 1

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

strActiveRMSE = cell(1,length(pVActiveRMSE));
for i = 1:length(pVActiveRMSE)
    if pVActiveRMSE(i)<0.001
        strActiveRMSE{i} = sprintf('%.3e', pVActiveRMSE(i));
    else
        strActiveRMSE{i} = sprintf('%.3f', pVActiveRMSE(i));
    end
end

RMSENN = readtable("...\ActiveTrig_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = median(RMSENN(RMSENN(:,1)~=0,Ind))';
load('RandomTrig.mat')
eOffline = median(RMSEOffline(RMSEOffline(:,1)~=0,Ind))';
fOffline = [quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(4))' ...
    quantile(RMSEOffline(RMSEOffline(:,1)~=0,Ind),qLevel(5))'];
load('ALTrig.mat')
eNT = median(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind))';
eST = median(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind))';
eLMC = median(RMSELMC(RMSELMC(:,1)~=0,Ind))';
eProp = median(RMSE(RMSE(:,1)~=0,Ind))';

fNN = [quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSENN(RMSENN(:,1)~=0,Ind),qLevel(5))')];
fNT = [quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind),qLevel(5))')];
fST = [quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind),qLevel(5))')];
fLMC = [quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(RMSELMC(RMSELMC(:,1)~=0,Ind),qLevel(5))')];
fProp = [quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(4))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(5))'];

figure
t = tiledlayout(4,2);
nexttile%([2 2]) % Mean RMSE of all methods for random batches
hold on
p1 = plot(eNT,...
    'LineStyle','-','Color',[255 0 255 255]/255,...
    'Marker','s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(eST,...
    'LineStyle','-','Color',[255 0 0 255]/255,...
    'Marker','d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(eLMC,...
    'LineStyle','-','Color',[197 90 17 255]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(eNN,...
    'LineStyle','-','Color',[153 67 96 255]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(eOffline,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0.2 1])
xticks([1 5:5:30])
yticks(0:0.2:1)
xticklabels('')
xlabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(a) Median RMSE for Active Batches','FontWeight','bold')
box on
grid on
% legend([p1,p2,p3,p4,p5,p6],{'SPGP','S2S','LMC','NN','Offline','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for SPGP
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNT(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNT(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('','FontSize',12)
title(['(b) SPGP (p=' strActiveRMSE{1} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for S2S
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fST(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fST(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('RMSE','FontSize',12)
title(['(c) S2S (p=' strActiveRMSE{2} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for LMC
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{i},'Color',[197 90 17 255]/255,...
            'LineWidth',1.5);
    else
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{end},'Color',[197 90 17 255]/255 ...
            ,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
ylabel('','FontSize',12)
title(['(d) LMC (p=' strActiveRMSE{3} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for NN
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNN(:,i),...
            'LineStyle',LineStyle{i},'Color',[153 67 96 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNN(:,i),...
            'LineStyle',LineStyle{end},'Color',[153 67 96 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
xlabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title(['(e) NN (p=' strActiveRMSE{4} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for Offline
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xticklabels('')
xlabel('','FontSize',12)
ylabel('','FontSize',12)
title(['(f) Offline (p=' strActiveRMSE{5} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of RMSE for Prop
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fProp(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fProp(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 2])
xticks([1 5:5:30])
yticks(0:0.5:2)
xlabel('Round','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(g) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

end


%% MAPE for Random Selection

for k = 1

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

strRandomMAPE = cell(1,length(pVActiveRMSE));
for i = 1:length(pVActiveRMSE)
    if pVRandomMAPE(i)<0.001
        strRandomMAPE{i} = sprintf('%.3e', pVRandomMAPE(i));
    else
        strRandomMAPE{i} = sprintf('%.3f', pVRandomMAPE(i));
    end
end

MAPENN = readtable("...\RandomTrig_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = median(MAPENN(MAPENN(:,1)~=0,Ind))';
load('RandomTrig.mat')
eOffline = median(MAPEOffline(MAPEOffline(:,1)~=0,Ind))';
% load('Random8S2D.mat')
eNT = median(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind))';
eST = median(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind))';
eLMC = median(MAPELMC(MAPELMC(:,1)~=0,Ind))';
eProp = median(MAPE(MAPE(:,1)~=0,Ind))';

fNN = [quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(5))')];
fNT = [quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(5))')];
fST = [quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(5))')];
fLMC = [quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(5))')];
fProp = [quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(4))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(5))'];
X = (1:25)';

figure
t = tiledlayout(4,2);
nexttile%([2 2]) % Median MAPE for all methods
hold on
p1 = plot(eNT,...
    'LineStyle','-','Color',[255 0 255 255]/255,...
    'Marker','s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(eST,...
    'LineStyle','-','Color',[255 0 0 255]/255,...
    'Marker','d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(eLMC,...
    'LineStyle','-','Color',[197 90 17 255]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(eNN,...
    'LineStyle','-','Color',[153 67 96 255]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(eOffline,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 150])
xticks([1 5:5:30])
yticks(0:25:150)
xticklabels('')
% xlabel('Round','FontSize',12)
ylabel('MAPE','FontSize',12)
title('(a) Median MAPE for Random Batches','FontWeight','bold')
box on
grid on
% legend({'NN','SPGP','S2S','LMC','Proposed'},'NumColumns',3)
% legend([p1,p2,p3,p4,p5],{'SPGP','S2S','LMC','NN','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

fOffline = [quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(5))')];
nexttile % Quantiles of MAPE for SPGP
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNT(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNT(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('','FontSize',12)
title(['(b) SPGP (p=' strRandomMAPE{1} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for S2S
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fST(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fST(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('MAPE','FontSize',12)
title(['(c) S2S (p=' strRandomMAPE{2} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for LMC
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{i},'Color',[197 90 17 255]/255,...
            'LineWidth',1.5);
    else
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{end},'Color',[197 90 17 255]/255 ...
            ,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('','FontSize',12)
title(['(d) LMC (p=' strRandomMAPE{3} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for NN
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNN(:,i),...
            'LineStyle',LineStyle{i},'Color',[153 67 96 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNN(:,i),...
            'LineStyle',LineStyle{end},'Color',[153 67 96 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 650])
xticks([1 5:5:30])
yticks(0:100:650)
xticklabels('')
xlabel('','FontSize',12)
ylabel('MAPE','FontSize',12)
title(['(e) NN (p=' strRandomMAPE{4} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for Offline
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('','FontSize',12)
ylabel('','FontSize',12)
title(['(f) Offline (p=' strRandomMAPE{5} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for Prop
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fProp(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fProp(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xlabel('Round','FontSize',12)
ylabel('MAPE','FontSize',12)
title('(g) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

end

%% MAPE for Active Selection

for k = 1

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

strActiveMAPE = cell(1,length(pVActiveRMSE));
for i = 1:length(pVActiveRMSE)
    if pVActiveMAPE(i)<0.001
        strActiveMAPE{i} = sprintf('%.3e', pVActiveMAPE(i));
    else
        strActiveMAPE{i} = sprintf('%.3f', pVActiveMAPE(i));
    end
end

MAPENN = readtable("...\ActiveTrig_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = median(MAPENN(MAPENN(:,1)~=0,Ind))';
load('RandomTrig.mat')
eOffline = median(MAPEOffline(MAPEOffline(:,1)~=0,Ind))';
load('ALTrig.mat')
eNT = median(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind))';
eST = median(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind))';
eLMC = median(MAPELMC(MAPELMC(:,1)~=0,Ind))';
eProp = median(MAPE(MAPE(:,1)~=0,Ind))';

fNN = [quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPENN(MAPENN(:,1)~=0,Ind),qLevel(5))')];
fNT = [quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind),qLevel(5))')];
fST = [quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind),qLevel(5))')];
fLMC = [quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPELMC(MAPELMC(:,1)~=0,Ind),qLevel(5))')];
fProp = [quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(4))' ...
    quantile(MAPE(MAPE(:,1)~=0,Ind),qLevel(5))'];
X = (1:25)';

figure
t = tiledlayout(4,2);
nexttile%([2 2]) % Median MAPE for all methods
hold on
p1 = plot(eNT,...
    'LineStyle','-','Color',[255 0 255 255]/255,...
    'Marker','s','LineWidth',1.5,'MarkerSize',MarkerSize);
p2 = plot(eST,...
    'LineStyle','-','Color',[255 0 0 255]/255,...
    'Marker','d','LineWidth',1.5,'MarkerSize',MarkerSize);
p3 = plot(eLMC,...
    'LineStyle','-','Color',[197 90 17 255]/255,...
    'Marker','<','LineWidth',1.5,'MarkerSize',MarkerSize);
p4 = plot(eNN,...
    'LineStyle','-','Color',[153 67 96 255]/255,...
    'Marker','^','LineWidth',1.5,'MarkerSize',MarkerSize);
p5 = plot(eOffline,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 100])
xticks([1 5:5:30])
yticks(0:20:100)
xticklabels('')
% xlabel('Round','FontSize',12)
ylabel('MAPE','FontSize',12)
title('(a) Median MAPE for Active Batches','FontWeight','bold')
box on
grid on
% legend({'NN','SPGP','S2S','LMC','Proposed'},'NumColumns',3)
% legend([p1,p2,p3,p4,p5],{'SPGP','S2S','LMC','NN','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

load('RandomTrig.mat')
fOffline = [quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(1))' ...
    quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(2))' ...
    quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(3))' ...
    quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(4))' ...
    (quantile(MAPEOffline(MAPEOffline(:,1)~=0,Ind),qLevel(5))')];
nexttile % Quantiles of MAPE for SPGP
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNT(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNT(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('','FontSize',12)
title(['(b) SPGP (p=' strActiveMAPE{1} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for S2S
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fST(:,i),...
            'LineStyle',LineStyle{i},'Color',[255 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fST(:,i),...
            'LineStyle',LineStyle{end},'Color',[255 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('MAPE','FontSize',12)
title(['(c) S2S (p=' strActiveMAPE{2} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for LMC
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{i},'Color',[197 90 17 255]/255,...
            'LineWidth',1.5);
    else
        plot(fLMC(:,i),...
            'LineStyle',LineStyle{end},'Color',[197 90 17 255]/255 ...
            ,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('','FontSize',12)
title(['(d) LMC (p=' strActiveMAPE{3} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for NN
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fNN(:,i),...
            'LineStyle',LineStyle{i},'Color',[153 67 96 255]/255,...
            'LineWidth',1.5);
    else
        plot(fNN(:,i),...
            'LineStyle',LineStyle{end},'Color',[153 67 96 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 600])
xticks([1 5:5:30])
yticks(0:100:600)
xticklabels('')
ylabel('MAPE','FontSize',12)
title(['(e) NN (p=' strActiveMAPE{4} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for Offline
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fOffline(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xticklabels('')
ylabel('','FontSize',12)
ylabel('','FontSize',12)
title(['(f) Offline (p=' strActiveMAPE{5} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for Prop
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fProp(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 255 255]/255,...
            'LineWidth',1.5);
    else
        plot(fProp(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 255 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0 250])
xticks([1 5:5:30])
yticks(0:50:250)
xlabel('Round','FontSize',12)
ylabel('MAPE','FontSize',12)
title('(g) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

end

