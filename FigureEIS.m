
clear
clc

%% Run code to generate results

GenerateResultsEISAL % This provides ALEIS.mat

% ALso please run the python code CaseEIS for NN

%% paired t-test

% Active batches
RMSENN = readtable("...\EIS_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = RMSENN(RMSENN(:,1)~=0,:);
eNN = eNN(:,end);
load('ALEIS.mat')
eNT = RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = RMSELMC(RMSELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = RMSE(RMSE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST eLMC eNN eProp];

pVActiveRMSE = zeros(1,4);
for i = 1:4
    [~,pVActiveRMSE(i)] = ttest(eAll(:,i),eAll(:,end));
end

% MAPE Active batches
MAPENN = readtable("...\EIS_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = MAPENN(MAPENN(:,1)~=0,:);
eNN = eNN(:,end);
load('ALEIS.mat')
eNT = MAPENonTransfer(MAPENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = MAPESGPTransfer(MAPESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = MAPELMC(MAPELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = MAPE(MAPE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST eLMC eNN eProp];

pVActiveMAPE = zeros(1,4);
for i = 1:4
    [~,pVActiveMAPE(i)] = ttest(eAll(:,i),eAll(:,end));
end

%% RMSE

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

strActiveRMSE = cell(1,length(pVActiveRMSE));
strActiveMAPE = cell(1,length(pVActiveRMSE));
for i = 1:length(pVActiveRMSE)
    if pVActiveRMSE(i)<0.001
        strActiveRMSE{i} = sprintf('%.3e', pVActiveRMSE(i));
    else
        strActiveRMSE{i} = sprintf('%.3f', pVActiveRMSE(i));
    end
    if pVActiveMAPE(i)<0.001
        strActiveMAPE{i} = sprintf('%.3e', pVActiveMAPE(i));
    else
        strActiveMAPE{i} = sprintf('%.3f', pVActiveMAPE(i));
    end
end

% MAPENN = readtable("D:\OneDrive - USTC\Research\Chao\05_Active learning\Code\NN4AL\Data4NN\EIS_MAPE.csv");
RMSENN = readtable("...\EIS_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = median(RMSENN(RMSENN(:,1)~=0,Ind))';
load('ALEIS.mat')
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
X = (1:25)';

figure
t = tiledlayout(3,2);
nexttile % Median MAPE for all methods
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
p5 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 0.1])
xticks([1 5:5:30])
yticks(0:0.02:0.1)
xticklabels('')
% xlabel('Round','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(a) Median RMSE','FontWeight','bold')
box on
grid on
% legend({'NN','SPGP','S2S','LMC','Proposed'},'NumColumns',3)
% legend([p1,p2,p3,p4,p5],{'SPGP','S2S','LMC','NN','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

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
axis([0.5 25.5 0 0.15])
xticks([1 5:5:30])
yticks(0:0.05:0.15)
xticklabels('')
title(['(b) SPGP (p=' strActiveRMSE{1} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 0.15])
xticks([1 5:5:30])
yticks(0:0.05:0.15)
xticklabels('')
ylabel('RMSE','FontSize',12)
title(['(c) S2S (p=' strActiveRMSE{2} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 0.15])
xticks([1 5:5:30])
yticks(0:0.05:0.15)
xticklabels('')
title(['(d) LMC (p=' strActiveRMSE{3} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 0.15])
xticks([1 5:5:30])
yticks(0:0.05:0.15)
% xticklabels('')
xlabel('Round','FontSize',12)
ylabel('RMSE','FontSize',12)
title(['(e) NN (p=' strActiveRMSE{4} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 0.15])
xticks([1 5:5:30])
yticks(0:0.05:0.15)
xlabel('Round','FontSize',12)
ylabel('')
title('(f) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';


%% MAPE comparison

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

strActiveRMSE = cell(1,length(pVActiveRMSE));
strActiveMAPE = cell(1,length(pVActiveRMSE));
for i = 1:length(pVActiveRMSE)
    if pVActiveRMSE(i)<0.001
        strActiveRMSE{i} = sprintf('%.3e', pVActiveRMSE(i));
    else
        strActiveRMSE{i} = sprintf('%.3f', pVActiveRMSE(i));
    end
    if pVActiveMAPE(i)<0.001
        strActiveMAPE{i} = sprintf('%.3e', pVActiveMAPE(i));
    else
        strActiveMAPE{i} = sprintf('%.3f', pVActiveMAPE(i));
    end
end

% MAPENN = readtable("D:\OneDrive - USTC\Research\Chao\05_Active learning\Code\NN4AL\Data4NN\EIS_MAPE.csv");
MAPENN = readtable("...\EIS_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = median(MAPENN(MAPENN(:,1)~=0,Ind))';
load('ALEIS.mat')
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
t = tiledlayout(3,2);
nexttile % Median MAPE for all methods
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
p5 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 40])
xticks([1 5:5:30])
yticks(0:10:40)
xticklabels('')
% xlabel('Round','FontSize',12)
ylabel('MAPE','FontSize',12)
title('(a) Median MAPE','FontWeight','bold')
box on
grid on
% legend({'NN','SPGP','S2S','LMC','Proposed'},'NumColumns',3)
% legend([p1,p2,p3,p4,p5],{'SPGP','S2S','LMC','NN','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

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
axis([0.5 25.5 0 80])
xticks([1 5:5:30])
yticks(0:20:80)
xticklabels('')
title(['(b) SPGP (p=' strActiveMAPE{1} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 80])
xticks([1 5:5:30])
yticks(0:20:80)
xticklabels('')
ylabel('MAPE','FontSize',12)
title(['(c) S2S (p=' strActiveMAPE{2} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 80])
xticks([1 5:5:30])
yticks(0:20:80)
xticklabels('')
title(['(d) LMC (p=' strActiveMAPE{3} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 80])
xticks([1 5:5:30])
yticks(0:20:80)
% xticklabels('')
xlabel('Round','FontSize',12)
ylabel('MAPE','FontSize',12)
title(['(e) NN (p=' strActiveMAPE{4} ')'],'FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for SPGP
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
axis([0.5 25.5 0 80])
xticks([1 5:5:30])
yticks(0:20:80)
xlabel('Round','FontSize',12)
ylabel('')
title('(f) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

