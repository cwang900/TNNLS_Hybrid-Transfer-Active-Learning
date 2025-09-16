
clear
clc

%% Run code to generate results

GenerateResultsSensitivity % This provides Sensitivity.mat

%%

for k = 1

MarkerSize = 4;
Ind = 2:26;
qLevel = 0.1:0.2:0.9;%[0.05 0.25 0.5 0.75 0.95];%
LineStyle = {':','-.','-','--','-','none'};
Marker = {'s','d','<','^','>','.'};

load('Sensitivity.mat')
eOL = mean(RMSEOnline(RMSEOnline(:,1)~=0,Ind))';
eOnline = median(RMSEOnline(RMSEOnline(:,1)~=0,Ind))';
fOnline = [quantile(RMSEOnline(RMSEOnline(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSEOnline(RMSEOnline(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSEOnline(RMSEOnline(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSEOnline(RMSEOnline(:,1)~=0,Ind),qLevel(4))' ...
    quantile(RMSEOnline(RMSEOnline(:,1)~=0,Ind),qLevel(5))'];

eP = mean(RMSE(RMSE(:,1)~=0,Ind))';
eProp = median(RMSE(RMSE(:,1)~=0,Ind))';
fProp = [quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(1))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(2))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(3))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(4))' ...
    quantile(RMSE(RMSE(:,1)~=0,Ind),qLevel(5))'];

figure
t = tiledlayout(3,1);
nexttile%([2 2]) % Mean RMSE of all methods for random batches
hold on
p5 = plot(eOL,...
    'LineStyle','-','Color',[0 0 0 255]/255,...
    'Marker','o','LineWidth',1.5,'MarkerSize',MarkerSize);
p6 = plot(eP,...
    'LineStyle','-','Color',[0 0 255 255]/255,...
    'Marker','>','LineWidth',1.5,'MarkerSize',MarkerSize);
hold off
axis([0.5 25.5 0.17 0.8])
xticks([1 5:5:30])
yticks(0.2:0.2:1)
xticklabels('')
xlabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(a) Mean RMSE','FontWeight','bold')
box on
grid on
legend([p6,p5],{'Proposed','Proposed with updating'},'NumColumns',3)
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

nexttile % Quantiles of MAPE for Online
hold on
for i = 1:length(qLevel)
    if i<length(qLevel)
        plot(fOnline(:,i),...
            'LineStyle',LineStyle{i},'Color',[0 0 0 255]/255,...
            'LineWidth',1.5);
    else
        plot(fOnline(:,i),...
            'LineStyle',LineStyle{end},'Color',[0 0 0 255]/255,...
            'Marker',Marker{end},'LineWidth',1.5,'MarkerSize',7,...
            "MarkerIndices",1:length(Ind));
    end
end
axis([0.5 25.5 0.17 1])
xticks([1 5:5:30])
yticks(0.2:0.2:1)
xticklabels('')
ylabel('','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(b) Proposed with Updating','FontWeight','bold')
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
axis([0.5 25.5 0.17 1])
xticks([1 5:5:30])
yticks(0.2:0.2:1)
xlabel('Round','FontSize',12)
ylabel('RMSE','FontSize',12)
title('(c) Proposed','FontWeight','bold')
box on
grid on
set(gca,"LineWidth",1.5,'FontWeight','bold','FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

end
