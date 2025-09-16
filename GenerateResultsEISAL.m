clear
clc

%%

NRep = 100;
LearningRound = 26;

NoiseSigma = [0.1 0.03];
RMSE = zeros(NRep,LearningRound);
RMSEOffline = zeros(NRep,LearningRound);
RMSESGPTransfer = zeros(NRep,LearningRound);
RMSENonTransfer = zeros(NRep,LearningRound);
RMSELMC = zeros(NRep,LearningRound);
MAPE = zeros(NRep,LearningRound);
MAPEOffline = zeros(NRep,LearningRound);
MAPESGPTransfer = zeros(NRep,LearningRound);
MAPENonTransfer = zeros(NRep,LearningRound);
MAPELMC = zeros(NRep,LearningRound);
tEnd = zeros(NRep,LearningRound);
tEndProp = zeros(NRep,LearningRound);
tEndSGPTransfer = zeros(NRep,LearningRound);
tEndNonTransfer = zeros(NRep,LearningRound);
tEndLMC = zeros(NRep,LearningRound);
FunVar = cell(1,5);

% parpool(15) % Set parallel computing
parfor iRep = 1:NRep
    %% Generate data

    NDomain = 3;
    NStream = 2;
    NObs = 30;
    y = cell(1,NDomain);
    X = cell(1,NDomain);

    %% Parametric simulation

    Data = load("NC_EIS_Capacity.mat");
    Data = Data.Data.T25;
    IndState = [1 2 5];
    IndCell = randperm(8,1);
    StateString = {'I','II','III','IV','V','VI','IX'};
    fAll = cell(1,8);
    XAll = cell(1,8);
    ii = 1;
    for iState = IndState
        Temp = ['Cell' num2str(IndCell)];
        State = ['State' StateString{iState}];
        fAll{ii} = [flipud(Data.(Temp).(State).EIS_Re(:,1)); ...
                flipud(Data.(Temp).(State).EIS_Im(:,1))];
        XAll{ii} = flipud(Data.(Temp).(State).Freq(:,1));
        XAll{ii} = repmat(XAll{ii},2,1);
        ii = ii+1;
    end
    X0 = XAll{1};
    Xmin = XAll{1}(1);
    Xmax = XAll{1}(end);

    for i = 1:NDomain
        if i<NDomain
            p = sort(lhsdesign(NObs,1));
            XX = Xmin+(Xmax-Xmin)*p;
            XTemp = reshape(XAll{i},[],2);
            yTemp = reshape(fAll{i},[],2);
            X{i} = [XX XX];
            X{i} = X{i}(:);
            ff = interp1(XTemp(:,1),yTemp,XX);
            ff = ff(:);
            y{i} = ff+normrnd(zeros(2*NObs,1),repelem(NoiseSigma',NObs,1));
        else
            p = sort(lhsdesign(5,1));
            XX = Xmin+(Xmax-Xmin)*p;
            XTemp = reshape(XAll{i},[],2);
            yTemp = reshape(fAll{i},[],2);
            X{i} = [XX XX];
            X{i} = X{i}(:);
            ff = interp1(XTemp(:,1),yTemp,XX);
            ff = ff(:);
            y{i} = ff+normrnd(zeros(10,1),repelem(NoiseSigma',5,1));
        end
    end
    ffun = [XAll{NDomain} fAll{NDomain}];

    %%

    XRand = unifrnd(Xmin,Xmax,2,LearningRound);
   
    CovType = {'CovOneShare','CovSGPTransfer','CovNonTransfer','CovLMC'};
    
    %% Set intial parameter
    
    NDomain = length(X);
    Rho0 = [1*unifrnd(1,1,1,NDomain-1);1*unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(3,3,1,NDomain-1);unifrnd(3,3,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1.5*[unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%/sqrt(NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([1*unifrnd(1,2,2,NDomain-1) 1*unifrnd(1,2,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%thetaSample;% [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%

    %% Active learning
    
    try
    [RMSE(iRep,:), MAPE(iRep,:), fTrue, XSample, ySample, tEndProp(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{1},LearningRound,XRand,NoiseSigma,NStream);
        
    NDomain = length(X);
    Rho0 = 1*[unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = 0.1*([unifrnd(5,5,1,NDomain-1);unifrnd(5,5,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*[unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%/sqrt(NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = 0.1*(1*[unifrnd(3,3,2,NDomain-1) 1*unifrnd(2,2,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 0.1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%thetaSample;% [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%

    [RMSESGPTransfer(iRep,:), MAPESGPTransfer(iRep,:), ~, ~, ~, tEndSGPTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{2},LearningRound,XRand,NoiseSigma,NStream);
    
    NDomain = length(X);
    Rho0 = 1*[unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = (1*[unifrnd(4,4,1,NDomain-1);unifrnd(4,4,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [unifrnd(1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([1*unifrnd(3,4,2,NDomain-1) 1*unifrnd(3,4,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(1,:)'; Scale(end); sigma(:)];
    [RMSEOffline(iRep,:), MAPEOffline(iRep,:), tEnd(iRep,:)] = OfflineLearning...
        (XSample,ySample,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,NStream);

    [RMSELMC(iRep,:), MAPELMC(iRep,:), ~, ~, ~, tEndLMC(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,XRand,NoiseSigma,NStream);
    
    NDomain = length(X);
    Rho0 = [unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(3,3,1,NDomain-1);unifrnd(3,3,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [unifrnd(1,1,2,NDomain-1) 0.1*unifrnd(1,1,2,1)];%/sqrt(NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(3,3,2,NDomain-1) 0.2*unifrnd(2,2,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 0.1*ones(1,NDomain);%

    y = y(end);
    X = X(end);
    Rho = Rho(:,end);%RhoSample(:,end);%unifrnd(1,2,2,NDomain);%
    Scale = Scale(:,end);%unifrnd(0.1,0.2,2,1);%ScaleSample(:,end);%
    sigma = sigma(1);
    theta0 = [Rho(:); Scale(:); sigma(:)];

    [RMSENonTransfer(iRep,:), MAPENonTransfer(iRep,:), ~, ~, ~, tEndNonTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{3},LearningRound,XRand,NoiseSigma,NStream);
    
    catch
        continue;
    end
end
save('ALEIS')

%%
function [RMSE, MAPE, fTrue, XSample, ySample, tEnd] = ActiveLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound,XRand,NoiseSigma,NStream)
tEnd = zeros(1,LearningRound);
NDomain = length(X);
% b = b(:,NDomain:NDomain:end);

%% Set input location for the integration

NObs = size(ffun,1)/2;
XStar = ffun(:,1);

NInd = 40;
XTInd = linspace(Xmin,Xmax,NInd)';
XTInd = [XTInd;XTInd];

%% Parameter optimization & Prediction

theta = OptJOffline(theta0,y,X,XTInd,CovType,NStream);

alpha = [];
C = [];
alphat = zeros(2*NInd,LearningRound);
Ct = zeros(2*NInd,2*NInd,LearningRound);
[Mu, Sigma, alphat(:,1), Ct(:,:,1)] = PostInd(theta,y,X,{XStar},alpha,C,{XTInd},NDomain,CovType,NStream);

%%

fStar1 = Mu(1:end/2);
fStar2 = Mu(end/2+1:end);

ff = ffun(:,2);
f1 = ff(1:end/2);
f2 = ff(end/2+1:end);
y1 = y{end}(1:end/2);
y2 = y{end}(end/2+1:end);
X1 = X{end}(1:end/2);
X2 = X{end}(end/2+1:end);
XStar1 = XStar(1:end/2);
XStar2 = XStar(end/2+1:end);
S = diag(Sigma);
S1 = S(1:end/2);
S2 = S(end/2+1:end);
% figure
% plot(X1,y1,'ko')
% hold on
% plot(X2,y2,'ko')
% plot(XStar1,f1,'b-x')
% plot(XStar2,f2,'b-+')
% plot(XStar1,fStar1,'r--')
% plot(XStar1,fStar1+2*sqrt(S1),'g--')
% plot(XStar1,fStar1-2*sqrt(S1),'g--')
% plot(XStar2,fStar2,'r--')
% plot(XStar2,fStar2+2*sqrt(S2),'g--')
% plot(XStar2,fStar2-2*sqrt(S2),'g--')

%% Active learning

XAL = X;
yAL = y;
NSample = 2;
iAL = 1;
XNew = [];
RMSE = zeros(1,LearningRound);
RMSE(1) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
RMSE(1) = sqrt(RMSE(1)/(2*NObs));

MAPE = zeros(1,LearningRound);
MAPE(1) = sum(abs((fStar1(:)-f1(:))./f1(:))+abs((fStar2(:)-f2(:))./f2(:)));
MAPE(1) = 100*(MAPE(1)/(2*NObs));

while iAL<LearningRound %Error > 1e-3
    tstart = tic;
    p = sort(lhsdesign(NSample,2),2);
    XSample0 = Xmin+p*(Xmax-Xmin);
    % XTemp = OPTIMSE(XSample0,XStar,Ct(:,:,iAL),XTInd,theta,NDomain,CovType,Xmin,Xmax,NStream);
    XTemp = XRand(:,iAL);
    XNew = [XNew XTemp(:)];

    ffNew = [interp1(XStar1,f1,XTemp) interp1(XStar2,f2,XTemp)];
    NAdd = size(ffNew,1);
    yNew = ffNew;%+normrnd(0,thetaAL(end,iAL),1,2);
    yNew = yNew+normrnd(zeros(size(yNew)),repelem(NoiseSigma,NAdd,1));
    yAL{end} = reshape(yAL{end},[],2);
    yAL{end} = [yAL{end}; yNew];
    yAL{end} = yAL{end}(:);
    XAL{end} = reshape(XAL{end},[],2);
    XAL{end} = [XAL{end}; [XTemp XTemp]];
    XAL{end} = XAL{end}(:);
    
    yt = yNew(:);
    Xt = [XTemp(:);XTemp(:)];

    tEnd(iAL) = toc(tstart);

    if strcmp(CovType,'CovNonTransfer')
        [Mu, Sigma, alphat(:,iAL+1), Ct(:,:,iAL+1)] = ...
            PostInd(theta,yAL,XAL,{XStar},alphat(:,iAL),Ct(:,:,iAL),{XTInd},NDomain,CovType,NStream);
    else
        [Mu, Sigma, alphat(:,iAL+1), Ct(:,:,iAL+1)] = ...
            PostInd(theta,{yt},{Xt},{XStar},alphat(:,iAL),Ct(:,:,iAL),{XTInd},NDomain,CovType,NStream);
    end
    iAL = iAL+1;

    fStar1 = Mu(1:end/2);
    fStar2 = Mu(end/2+1:end);
    y1 = yAL{end}(1:end/2);
    y2 = yAL{end}(end/2+1:end);
    X1 = XAL{end}(1:end/2);
    X2 = XAL{end}(end/2+1:end);
    S = diag(Sigma);
    S1 = S(1:end/2);
    S2 = S(end/2+1:end);
    % figure
    % plot(X1,y1,'ko')
    % hold on
    % plot(X2,y2,'ko')
    % plot(XNew,0*XNew,'kx','MarkerSize',12)
    % plot(XStar1,f1,'b-x')
    % plot(XStar2,f2,'b-+')
    % plot(XStar1,fStar1,'r--')
    % plot(XStar1,fStar1+2*sqrt(S1),'g--')
    % plot(XStar1,fStar1-2*sqrt(S1),'g--')
    % plot(XStar2,fStar2,'r--')
    % plot(XStar2,fStar2+2*sqrt(S2),'g--')
    % plot(XStar2,fStar2-2*sqrt(S2),'g--')

    RMSE(iAL) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
    RMSE(iAL) = sqrt(RMSE(iAL)/(2*NObs));

    MAPE(iAL) = sum(abs((fStar1(:)-f1(:))./f1(:))+abs((fStar2(:)-f2(:))./f2(:)));
    MAPE(iAL) = 100*(MAPE(iAL)/(2*NObs));
    
end

if strcmp(CovType,'CovOneShare')
    fTrue = [f1 f2];
    XSample = XAL;
    ySample = yAL;
else
    fTrue = [];
    XSample = [];
    ySample = [];
end

end

%%

function [Error, tEnd] = OfflineLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound)

%% Set input location for the integration

NObs = size(ffun,1)/2;
XStar = ffun(:,1);

ff = ffun(:,2);
f1 = ff(1:end/2);
f2 = ff(end/2+1:end);

Error = zeros(1,LearningRound);
tEnd = zeros(1,LearningRound);

for i = 0:LearningRound-1
    
    NDomain = length(y);
    tStart = tic;
    yTemp = y;
    yTemp{end} = reshape(yTemp{end},[],2);
    yTemp{end}(end-2*(LearningRound-1-i)+1:end,:) = [];
    yTemp{end} = yTemp{end}(:);
    XTemp = X;
    XTemp{end} = reshape(XTemp{end},[],2);
    XTemp{end}(end-2*(LearningRound-1-i)+1:end,:) = [];
    XTemp{end} = XTemp{end}(:);
    theta = OptJOffline(theta0,yTemp,XTemp,[],CovType);
    [~, ~, Mu, Sigma] = PostInd(theta,yTemp,XTemp,[],[],[],{XStar},NDomain,CovType);
    % [Mu, Sigma] = PostOffline(theta,yTemp,XTemp,XStar,[],CovType);

    fStar1 = Mu(1:end/2);
    fStar2 = Mu(end/2+1:end);

    Error(i+1) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
    Error(i+1) = sqrt(Error(i+1)/(2*NObs));
    
    tEnd(i+1) = toc(tStart);
end

end