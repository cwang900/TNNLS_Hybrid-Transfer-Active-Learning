clear
clc

%% Run code to generate results

GenerateResultsGrapheneAL % This provides ALGraphene.mat

% ALso please run the python code CaseGraphene for NN

%% Demonstration of the curves

Data = load("Graphene4Fig.mat");
Data = Data.Graphene4Fig;
fAll = cell(1,3);
XAll = cell(1,3);
NObs = 30;
NoiseSigma = [0.1 0.1];
Color = {[84 130 53]/255,'r','b'};
ii = 1;
figure
hold on
for iState = [1 3 2]%1:3
    if ii < 3
        NObs = 30;
    else
        NObs = 5;
    end
    fAll{ii} = [Data.f{1}((iState-1)*2+1,:)'; ...
        Data.f{1}((iState-1)*2+2,:)'];
    XAll{ii} = repmat(Data.X{1}',2,1);
    plot(XAll{ii}(1:end/2),fAll{ii}(1:end/2),...
        'Color',Color{ii},"LineWidth",1)
    plot(XAll{ii}(end/2+1:end),fAll{ii}(end/2+1:end),...
        'Color',Color{ii},'LineStyle','--',"LineWidth",1)

    Xmin = XAll{1}(1);
    Xmax = XAll{1}(end);

    p = sort(lhsdesign(NObs,1));
    XX = Xmin+(Xmax-Xmin)*p;
    XTemp = reshape(XAll{ii},[],2);
    yTemp = reshape(fAll{ii},[],2);
    X = [XX XX];
    X = X(:);
    ff = interp1(XTemp(:,1),yTemp,XX);
    ff = ff(:);
    y = ff+normrnd(zeros(2*NObs,1),repelem(NoiseSigma',NObs,1));
    plot(X(1:end/2),y(1:end/2),'Color',Color{ii},...
        'LineStyle','none','Marker','x','MarkerSize',7,"LineWidth",1.5)
    plot(X(end/2+1:end),y(end/2+1:end),'Color',Color{ii},...
        'LineStyle','none','Marker','o','MarkerSize',7,"LineWidth",1.5)

    ii = ii+1;
end
axis([0.5 2.5 2 20])
xticks(0.5:0.5:2.5)
yticks(0:5:30)
xlabel(['\it V' '\rm _{gs}\bf(V)'],'FontSize',12)
ylabel(['\it I' '\rm _{ds}\bf(\muA)'],'FontSize',12)
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)


%% Mean RMSE

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

figure
t = tiledlayout(3,2);
RMSENN = readtable("...\Graphene_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = mean(RMSENN(RMSENN(:,1)~=0,Ind))';
load('ALGraphene.mat')
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
p6 = plot(eProp,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0 1.5])
xticks([1 5:5:30])
yticks(0:0.3:1.5)
xticklabels('')
xlabel('','FontSize',12,'FontWeight','bold')
ylabel('RMSE','FontSize',12,'FontWeight','bold')
title('(a) RMSE for AL batches','FontWeight','bold')
box on
grid on
legend([p1,p2,p3,p4,p6],{'SPGP','S2S','LMC','NN','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

%% Plot of computational time for all online methods

tEndNN = readtable("...\ActiveGraphene_Time.csv");
tEndNN = tEndNN{3:27,:};
tEndNN = mean(tEndNN(tEndNN(:,2)~=0,:),2);
load('ALGraphene.mat')
tEndNonTransfer= mean(tEndNonTransfer(tEndNonTransfer(:,2)~=0,Ind-1));
tEndSGPTransfer = mean(tEndSGPTransfer(tEndSGPTransfer(:,2)~=0,Ind-1));
tEndLMC = mean(tEndLMC(tEndLMC(:,2)~=0,Ind-1));
tEndProp = mean(tEndProp(tEndProp(:,2)~=0,Ind-1));

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
axis([0.5 25.5 0 0.8])
xticks([1 5:5:30])
yticks([0:0.2:0.8])
% xticklabels('')
xlabel('Round','FontSize',12)
ylabel('Time (s)','FontSize',12)
box on
grid on
title('(b) Computational time','FontSize',12,'FontWeight','bold')
% legend([p1,p2,p3,p4,p5,p6],{'SPGP','S2S','LMC','NN','Offline','Proposed'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)%,'Offline'

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% Paired t-test

% Active batches
RMSENN = readtable("C:\Users\hzyll\OneDrive - USTC\Research\Chao\05_Active learning\Code\NN4AL\Data4NN\Graphene_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = RMSENN(RMSENN(:,1)~=0,:);
eNN = eNN(:,end);
load('ALResultsGrapheneOffline.mat')
eNT = RMSENonTransfer(RMSENonTransfer(:,1)~=0,Ind);
eNT = eNT(:,end);
eST = RMSESGPTransfer(RMSESGPTransfer(:,1)~=0,Ind);
eST = eST(:,end);
eLMC = RMSELMC(RMSELMC(:,1)~=0,Ind);
eLMC = eLMC(:,end);
eProp = RMSE(RMSE(:,1)~=0,Ind);
eProp = eProp(:,end);
eAll = [eNT eST(1:98) eLMC eNN(1:98) eProp(1:98)];

pVActiveRMSE = zeros(1,4);
for i = 1:4
    [~,pVActiveRMSE(i)] = ttest(eAll(:,i),eAll(:,end));
end

% MAPE Active batches
MAPENN = readtable("C:\Users\hzyll\OneDrive - USTC\Research\Chao\05_Active learning\Code\NN4AL\Data4NN\Graphene_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = MAPENN(MAPENN(:,1)~=0,:);
eNN = eNN(:,end);
load('ALResultsGraphene.mat')
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

%% Quantiles of RMSE

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

RMSENN = readtable("...\Graphene_RMSE.csv");
RMSENN = RMSENN{2:27,:}';
eNN = median(RMSENN(RMSENN(:,1)~=0,Ind))';
load('ALGraphene.mat')
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
axis([0.5 25.5 0 0.8])
xticks([1 5:5:30])
yticks(0:0.2:0.8)
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
axis([0.5 25.5 0 3])
xticks([1 5:5:30])
yticks(0:0.5:3)
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
axis([0.5 25.5 0 3])
xticks([1 5:5:30])
yticks(0:0.5:3)
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
axis([0.5 25.5 0 3])
xticks([1 5:5:30])
yticks(0:0.5:3)
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
axis([0.5 25.5 0 3])
xticks([1 5:5:30])
yticks(0:0.5:3)
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
axis([0.5 25.5 0 3])
xticks([1 5:5:30])
yticks(0:0.5:3)
xlabel('Round','FontSize',12)
ylabel('')
title('(f) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% MAPE

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

MAPENN = readtable("...\Graphene_MAPE.csv");
MAPENN = MAPENN{2:27,:}';
eNN = median(MAPENN(MAPENN(:,1)~=0,Ind))';
load('ALGraphene.mat')
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
axis([0.5 25.5 0 5])
xticks([1 5:5:30])
yticks(0:1:5)
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
axis([0.5 25.5 0 20])
xticks([1 5:5:30])
yticks(0:5:20)
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
axis([0.5 25.5 0 20])
xticks([1 5:5:30])
yticks(0:5:20)
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
axis([0.5 25.5 0 20])
xticks([1 5:5:30])
yticks(0:5:20)
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
axis([0.5 25.5 0 20])
xticks([1 5:5:30])
yticks(0:5:20)
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
axis([0.5 25.5 0 20])
xticks([1 5:5:30])
yticks(0:5:20)
xlabel('Round','FontSize',12)
ylabel('')
title('(f) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

