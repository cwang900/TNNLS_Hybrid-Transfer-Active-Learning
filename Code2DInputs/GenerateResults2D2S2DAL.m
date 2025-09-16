clear
clc
% close all

%%

NRep = 100;
LearningRound = 26;
NoiseSigma = 0.1;
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

% parpool(10) % Set parallel computing
parfor iRep = 1:NRep
    %% Generate data

    NDomain = 2;
    NStream = 2;
    NObs = 15;
    NObsT = 7;
    Xmin = 0;
    Xmax = 5;
    X0 = linspace(Xmin,Xmax,1000)';
    y = cell(1,NDomain);
    X = cell(1,NDomain);

    %% Parametric simulation
    
    for i = 1:NDomain
        if i<NDomain
            p = sort(lhsdesign(NObs,2));
            [x1, x2] = meshgrid(p(:,1),p(:,2));
            X{i} = Xmin+(Xmax-Xmin)*[x1(:) x2(:)];
        else
            p = linspace(Xmin,Xmax,NObsT)';%Xmin+(Xmax-Xmin)*sort(lhsdesign(NObsT,NStream));%
            [x1, x2] = meshgrid(p(:,1),p(:,1));
            X{i} = [x1(:) x2(:)];
        end
    end

    w0 = 1.5;
    rho = [[0.6; 0.4] [0.4; 0.6] [0.3; 0.7]];
    g1 = @(u,x,rho) u+3*((1-rho)*sin(w0*x(:,1))+rho*cos(w0*x(:,2))+cos(0.5*x(:,2)));
    g2 = @(u,x,rho) u+3*((1-rho)*cos(w0*x(:,2))+rho*sin(w0*x(:,1))+sin(0.5*x(:,1)));
    g = @(u,x,eta,rho) u+...
        exp(-(x(:,1)+x(:,2))./eta).*(rho*2*sin(w0*x(:,1))+(1-rho)*2*cos(w0*x(:,2)));

    Fun = {g1,g2};
    eta = unifrnd(-10,-5,NStream,1);%[unifrnd(-10,-5,NStream/2,1);unifrnd(4,5,NStream/2,1)];%[5,-10];
    U = unifrnd(5,10,NStream,NDomain);%zeros(2,NDomain);%
    
    f = cell(1,NDomain);
    for i = 1:NDomain
        f{i} = [];
        if i<NDomain
            f{i} = [Fun{i}(U(1,i),X{i},rho(1,i)); Fun{i}(U(2,i),X{i},rho(2,i))];
        else
            for ii = 1:NStream
                f{i} = [f{i}; g(U(ii,i),X{i},eta(ii),rho(ii,i))];
            end
        end
        X{i} = repmat(X{i},NStream,1);
        y{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
    end
    ffun = @(t)TrueFun(t,U(:,NDomain),eta,rho(:,NDomain),g);

    %%

    XRand = unifrnd(Xmin,Xmax,NStream,100);
   
    CovType = {'CovCP2D','CovSGPTransfer2D','CovNonTransfer2D','CovLMC2D'};
    
    %% 

    try

    %% Proposed
    NDomain = length(X);
    Rho0 = 1*unifrnd(1,1,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = 1*[3;1;1;3];%0.1*[[3;3;3;3] [3;3;3;3]];%unifrnd(1,1,2*NStream,NDomain-1);
    Rho = 1*[unifrnd(1,1,NStream,NDomain-1) 1*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = [1.5*[1.5;2;2;1.5] 0.2*[0.5;0.5;0.5;0.5]];%0.1*[[4;4;4;4] [4;4;4;4] [4;4;4;4]];%0.5*[[3;4;3;4] [4;3;4;3] [3;3;3;3]];%
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%thetaSample;% 

    NStream = repmat(NStream,NDomain,1);
    [RMSE(iRep,:), MAPE(iRep,:), fTrue, XSample, ySample, tEndProp(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{1},LearningRound,XRand,NoiseSigma,NStream);

    %% S2S

    NStream = NStream(1);
    NDomain = length(X);
    Rho0 = 1*unifrnd(1,1,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = 0.2*[3;4;3;4];%0.19*[[3;4;3;4] [4;3;4;3]];%unifrnd(1,1,2*NStream,NDomain-1);
    Rho = [1*unifrnd(1,1,NStream,NDomain-1) 1*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = 1.3*[[3;4;3;4] [4;3;4;3]];%0.15*[[3;4;3;4] [4;3;4;3] [3;3;3;3]];%([unifrnd(4,4,2*NStream,NDomain-1) unifrnd(2,2,2*NStream,1);]);
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];
    NStream = repmat(NStream,NDomain,1);
    [RMSESGPTransfer(iRep,:), MAPESGPTransfer(iRep,:), ~, ~, ~, tEndSGPTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{2},LearningRound,XRand,NoiseSigma,NStream);

    %% LMC

    NStream = NStream(1);
    NDomain = length(X);
    NLatent = 1;
    Rho0 = 1*unifrnd(1,1,NStream*NLatent,NDomain-1);
    Scale0 = 0.2*[[5;5] [5;5]];% 0.5*unifrnd(1,1,2,NLatent*NDomain);
    Rho = [1*unifrnd(1,1,NStream*NLatent,NDomain-1) 1*unifrnd(1,1,NStream*NLatent,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = 4*([unifrnd(0.5,0.5,NStream,NDomain-1) unifrnd(0.5,0.5,NStream,1)]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); sigma(:)];
    NStream = repmat(NStream,NDomain,1);
    [RMSEOffline(iRep,:), MAPEOffline(iRep,:), tEnd(iRep,:)] = OfflineLearning...
        (XSample,ySample,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,NStream);

    NStream = NStream(1);
    NDomain = length(X);
    NLatent = 1;
    Rho0 = 1*unifrnd(1,1,NStream*NLatent,NDomain-1);
    Scale0 = 0.2*[[5;5] [5;5]];% 0.5*unifrnd(1,1,2,NLatent*NDomain);
    Rho = [1*unifrnd(1,1,NStream*NLatent,NDomain-1) 1*unifrnd(1,1,NStream*NLatent,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = 4*([unifrnd(0.5,0.5,NStream,NDomain-1) unifrnd(0.5,0.5,NStream,1)]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 0.1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0(:); sigma(:)];
    NStream = repmat(NStream,NDomain,1);

    [RMSELMC(iRep,:), MAPELMC(iRep,:), ~, ~, ~, tEndLMC(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,XRand,NoiseSigma,NStream);

    %% SPGP

    NStream = NStream(1);
    Rho0 = 1*unifrnd(1,1,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = [[3;4;3;4] [4;3;4;3]];%unifrnd(1,1,2*NStream,NDomain-1);
    Rho = [1*unifrnd(1,1,NStream,NDomain-1) 1*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = 0.3*[[3;4;3;4] [4;3;4;3] [3;3;3;3]];%([unifrnd(4,4,2*NStream,NDomain-1) unifrnd(2,2,2*NStream,1);]);
    sigma = 1*ones(1,NDomain);%

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
save('AL2D2S2D.mat')

%%
function f = TrueFun(t,u,eta,rho,g)
    NStream = length(u);
    N = size(t,1)/NStream;
    f = [];
    for i = 1:NStream
        f = [f g(u(i),t((i-1)*N+1:i*N,:),eta(i),rho(i))];
    end
end

%%
function [RMSE, MAPE, fTrue, XSample, ySample, tEnd] = ActiveLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound,XRand,NoiseSigma,NStream)
tEnd = zeros(1,LearningRound);
NDomain = length(X);
D = size(X{1},2);

%% Set input location for the integration

NObs = 20;
tt = linspace(Xmin,Xmax,NObs)';
[X1, X2] = meshgrid(tt,tt);
XStar = repmat([X1(:) X2(:)],NStream(NDomain),1);

NInd = 10;
XTInd = linspace(Xmin,Xmax,NInd)';
[X1, X2] = meshgrid(XTInd,XTInd);
XTInd = repmat([X1(:) X2(:)],NStream(NDomain),1);

%% Parameter optimization & Prediction
% tic
theta = OptJOffline2D(theta0,y,X,XTInd,CovType,NStream);
% toc

alpha = [];
C = [];
% alphat = zeros(NStream(end)*NInd^2,LearningRound);
% Ct = zeros(NStream(end)*NInd^2,NStream(end)*NInd^2,LearningRound);
[Mu, Sigma, alphat, Ct] = PostInd2D(theta,y,X,{XStar},alpha,C,{XTInd},NDomain,CovType,NStream);

%%

fStar = Mu;
ff = ffun(XStar);

%% Active learning

XAL = X;
yAL = y;
NSample = 2;
iAL = 1;
XNew = [];
RMSE = zeros(1,LearningRound);
RMSE(1) = sum((fStar(:)-ff(:)).^2);
RMSE(1) = sqrt(RMSE(1)/(2*NObs^2));

MAPE = zeros(1,LearningRound);
MAPE(1) = sum(abs((fStar(:)-ff(:))./ff(:)));
MAPE(1) = 100*(MAPE(1)/(2*NObs^2));
% figure(1)
% hold on
% surf(reshape(XStar(1:NObs^2,1),NObs,NObs),reshape(XStar(1:NObs^2,2),NObs,NObs),...
%     reshape(ff(:,1),NObs,NObs))
% plot3(X{end}(1:end/2,1),X{end}(1:end/2,2),y{end}(1:end/2,1),'ko')
% surf(reshape(XStar(NObs^2+1:end,1),NObs,NObs),reshape(XStar(NObs^2+1:end,2),NObs,NObs),...
% reshape(fStar(1:NObs^2),NObs,NObs))
% grid on
% figure(2)
% hold on
% surf(reshape(XStar(NObs^2+1:end,1),NObs,NObs),reshape(XStar(NObs^2+1:end,2),NObs,NObs),...
%     reshape(ff(:,2),NObs,NObs))
% plot3(X{end}(1:end/2,1),X{end}(1:end/2,2),y{end}(end/2+1:end,1),'ko')
% surf(reshape(XStar(NObs^2+1:end,1),NObs,NObs),reshape(XStar(NObs^2+1:end,2),NObs,NObs),...
% reshape(fStar(NObs^2+(1:400)),NObs,NObs))
% grid on
while iAL<LearningRound %Error > 1e-3
    tstart = tic;
    p = sort(lhsdesign(NSample,D),1);
    XSample0 = unifrnd(Xmin,Xmax,2,NSample);%Xmin+p'*(Xmax-Xmin);
    XTemp = OPTIMSE2D(XSample0,XStar,Ct,XTInd,theta,NDomain,CovType,Xmin,Xmax,NStream);
    XNew = [XNew XTemp(:)];
    XTemp = reshape(XTemp(:)',[],D);
    ffNew = ffun(repmat(XTemp,NStream(end),1));
    % figure(1)
    % plot3(XTemp(:,1),XTemp(:,2),ffNew(:,1),'ro','MarkerFaceColor','r')
    % plot3(XTemp(:,1),XTemp(:,2),0*ffNew(:,1),'kx')
    % figure(2)
    % plot3(XTemp(:,1),XTemp(:,2),ffNew(:,2),'ro','MarkerFaceColor','r')
    % plot3(XTemp(:,1),XTemp(:,2),0*ffNew(:,2),'kx')
    yNew = ffNew;%+normrnd(0,thetaAL(end,iAL),1,2);
    yNew = yNew+normrnd(0,NoiseSigma,size(yNew));
    yAL{end} = reshape(yAL{end},[],NStream(end));
    yAL{end} = [yAL{end}; yNew];
    yAL{end} = yAL{end}(:);
    XAL{end} = XAL{end}(1:end/2,:);
    XAL{end} = [XAL{end}; XTemp];
    XAL{end} = repmat(XAL{end},NStream(end),1);
    
    yt = yNew(:);
    Xt = repmat(XTemp,NStream(end),1);
    
    tEnd(iAL) = toc(tstart);

    if strcmp(CovType,'CovNonTransfer')
        [Mu, Sigma, alphat, Ct] = ...
            PostInd2D(theta,yAL,XAL,{XStar},alphat,Ct,{XTInd},NDomain,CovType,NStream);
    else
        [Mu, Sigma, alphat, Ct] = ...
            PostInd2D(theta,{yt},{Xt},{XStar},alphat,Ct,{XTInd},NDomain,CovType,NStream);
    end
    iAL = iAL+1;
    
    fStar = Mu;

    RMSE(iAL) = sum((fStar(:)-ff(:)).^2);
    RMSE(iAL) = sqrt(RMSE(iAL)/(2*NObs^2));

    MAPE(iAL) = sum(abs((fStar(:)-ff(:))./ff(:)));
    MAPE(iAL) = 100*(MAPE(iAL)/(2*NObs^2));

end

if strcmp(CovType,'CovCP2D')
    fTrue = ff;
    XSample = XAL;
    ySample = yAL;
else
    fTrue = [];
    XSample = [];
    ySample = [];
end

end
