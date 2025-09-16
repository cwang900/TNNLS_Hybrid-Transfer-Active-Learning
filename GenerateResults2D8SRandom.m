clear
clc
% close all

%%

NRep = 100;
LearningRound = 30;
NoiseSigma = 0.5;
Error = zeros(NRep,LearningRound);
ErrorOffline = zeros(NRep,LearningRound);
ErrorSGPTransfer = zeros(NRep,LearningRound);
ErrorNonTransfer = zeros(NRep,LearningRound);
ErrorLMC = zeros(NRep,LearningRound);
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
    NStream = 8;
    NObs = 25;
    Xmin = -5;
    Xmax = 5;
    X0 = linspace(Xmin,Xmax,1000)';
    y = cell(1,NDomain);
    X = cell(1,NDomain);

    %% Parametric simulation
    
    NSample = NObs;
    p = sort(lhsdesign(NSample,NDomain));
    t = Xmin+(Xmax-Xmin)*p;
    t = mat2cell(t,NSample,ones(1,NDomain));
    % rng(1234)
    p = sort(lhsdesign(5,1));
    % rng("shuffle")
    t{end} = Xmin+(Xmax-Xmin)*p;

    w0 = 0.7;
    rho = lhsdesign(NStream,NDomain);
    g1 = @(u,x,rho) u+5*((1-rho)*sin(w0*x)+rho*sin(1*x));
    g2 = @(u,x,rho) u+5*((1-rho)*cos(w0*x)+rho*sin(1*x));
    g = @(u,x,eta,rho) u+...
        exp(-x./eta).*(rho*2*sin(w0*x)+(1-rho)*2*cos(w0*x));

    Fun = {g1,g2};
    eta = [unifrnd(-10,-5,NStream,1);unifrnd(4,5,NStream,1)];%[5,-10];
    U = unifrnd(-5,5,NStream,NDomain);%zeros(2,NDomain);%
    
    f = cell(1,NDomain);
    for i = 1:NDomain
        X{i} = repmat(t{i},NStream,1);
        f{i} = [];
        if i<NDomain
            for ii = 1:NStream
                f{i} = [f{i}; g1(U(ii,i),t{i},rho(ii,i))];
            end
        else
            for ii = 1:NStream
                f{i} = [f{i}; g(U(ii,i),t{i},eta(ii),rho(ii,i))];
            end
        end
        y{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
    end
    ffun = @(t)TrueFun(t,U(:,NDomain),eta,rho(:,NDomain),g);

    %%

    XRand = unifrnd(Xmin,Xmax,2,LearningRound);
   
    CovType = {'CovOneShare','CovSGPTransfer','CovNonTransfer','CovLMC'};
    % FunVar{5} = CovType;

    %% Set intial parameter

    NDomain = length(X);
    Rho0 = 1*unifrnd(1,1,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = unifrnd(1,1,NStream,NDomain-1);%unifrnd(2,3,2,NDomain-1);%
    Rho = [1*unifrnd(1,1,NStream,NDomain-1) 1*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(4,4,NStream,NDomain-1) unifrnd(2,2,NStream,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%thetaSample;% 

    %% Active learning
    
    try
    [Error(iRep,:), fTrue, XSample, ySample, tEndProp(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{1},LearningRound,XRand,NoiseSigma,NStream);

    NDomain = length(X);
    Rho0 = 1*unifrnd(1,1,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = unifrnd(1,1,NStream,NDomain-1);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*[unifrnd(1,1,NStream,NDomain-1) 1*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(5,5,NStream,NDomain-1) unifrnd(5,5,NStream,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];
    % 
    [ErrorSGPTransfer(iRep,:), ~, ~, ~, tEndSGPTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{2},LearningRound,XRand,NoiseSigma,NStream);

    NDomain = length(X);
    Rho0 = 1*unifrnd(1,1.5,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = 1*unifrnd(1,1,1,NDomain-1);%unifrnd(2,3,2,NDomain-1);%
    Rho = [1*unifrnd(1,1,NStream,NDomain) 0.5*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(1,2,1,1) unifrnd(0.1,0.2,1,1)]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0'; Scale'; sigma(:)];
    [ErrorOffline(iRep,:), tEnd(iRep,:)] = OfflineLearning...
        (XSample,ySample,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,NStream);
    
    NDomain = length(X);
    Rho0 = 1*unifrnd(1,1.5,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = 1*unifrnd(1,1,1,NDomain-1);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*unifrnd(1,1,NStream,NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = unifrnd(1,2,1,1);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    sigma = 1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0'; Scale(end)'; sigma(:)];

    [ErrorLMC(iRep,:), ~, ~, ~, tEndLMC(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,XRand,NoiseSigma,NStream);

    NDomain = length(X);
    Rho0 = unifrnd(1,1,NStream,NDomain-1);%unifrnd(1,3,2,NDomain-1);%
    Scale0 = unifrnd(4,4,NStream,NDomain-1);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*[unifrnd(1,1,NStream,NDomain-1) 1*unifrnd(1,1,NStream,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(1,1,NStream,NDomain-1) unifrnd(1,1,NStream,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 1*ones(1,NDomain);%

    y = y(end);
    X = X(end);
    Rho = Rho(:,end);%RhoSample(:,end);%unifrnd(1,2,2,NDomain);%
    Scale = Scale(:,end);%unifrnd(0.1,0.2,2,1);%ScaleSample(:,end);%
    sigma = sigma(1);
    theta0 = [Rho(:); Scale(:); sigma(:)];
    [ErrorNonTransfer(iRep,:), ~, ~, ~, tEndNonTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{3},LearningRound,XRand,NoiseSigma,NStream);

    catch
        continue;
    end
end
save('Random2D8S.mat')

%%

function f = TrueFun(t,u,eta,rho,g)
    NStream = length(u);
    f = [];
    for i = 1:NStream
        f = [f g(u(i),t,eta(i),rho(i))];
    end
end

%%
function [Error, fTrue, XSample, ySample, tEnd] = ActiveLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound,XRand,NoiseSigma,NStream)
tEnd = zeros(1,LearningRound);
NDomain = length(X);
% b = b(:,NDomain:NDomain:end);

%% Set input location for the integration

NObs = 50;
tt = linspace(Xmin,Xmax,NObs)';
XStar = repmat(tt,NStream,1);

NInd = 40;
XTInd = linspace(Xmin,Xmax,NInd)';
XTInd = repmat(XTInd,NStream,1);

%% Parameter optimization & Prediction

theta = OptJOffline(theta0,y,X,XTInd,CovType,NStream);

alpha = [];
C = [];
alphat = zeros(NStream*NInd,LearningRound);
Ct = zeros(NStream*NInd,NStream*NInd,LearningRound);
[Mu, Sigma, alphat(:,1), Ct(:,:,1)] = PostInd(theta,y,X,{XStar},alpha,C,{XTInd},NDomain,CovType,NStream);

%%

fStar = Mu;
ff = ffun(tt);

%% Active learning

XAL = X;
yAL = y;
NSample = 2;
iAL = 1;
XNew = [];
Error = zeros(1,LearningRound);
Error(1) = sum((fStar(:)-ff(:)).^2);
Error(1) = sqrt(Error(1)/(2*NObs));

while iAL<LearningRound %Error > 1e-3
    tstart = tic;
    p = sort(lhsdesign(NSample,2),2);
    XSample0 = Xmin+p*(Xmax-Xmin);
    % XTemp = OPTIMSE(XSample0,XStar,Ct(:,:,iAL),XTInd,theta,NDomain,CovType,Xmin,Xmax,NStream);
    XTemp = XRand(:,iAL);
    XNew = [XNew XTemp(:)];
    ffNew = ffun(XTemp(:));
    yNew = ffNew;%+normrnd(0,thetaAL(end,iAL),1,2);
    yNew = yNew+normrnd(0,NoiseSigma,size(yNew));
    yAL{end} = reshape(yAL{end},[],NStream);
    yAL{end} = [yAL{end}; yNew];
    yAL{end} = yAL{end}(:);
    XAL{end} = reshape(XAL{end},[],NStream);
    XAL{end} = [XAL{end}; repmat(XTemp,1,NStream)];
    XAL{end} = XAL{end}(:);
    
    yt = yNew(:);
    Xt = repmat(XTemp(:),NStream,1);
    
    tEnd(iAL) = toc(tstart);

    if strcmp(CovType,'CovNonTransfer')
        [Mu, Sigma, alphat(:,iAL+1), Ct(:,:,iAL+1)] = ...
            PostInd(theta,yAL,XAL,{XStar},alphat(:,iAL),Ct(:,:,iAL),{XTInd},NDomain,CovType,NStream);
    else
        [Mu, Sigma, alphat(:,iAL+1), Ct(:,:,iAL+1)] = ...
            PostInd(theta,{yt},{Xt},{XStar},alphat(:,iAL),Ct(:,:,iAL),{XTInd},NDomain,CovType,NStream);
    end
    iAL = iAL+1;
    
    fStar = Mu;
    % fStar1 = Mu(1:end/2);
    % fStar2 = Mu(end/2+1:end);
    % y1 = yAL{end}(1:end/2);
    % y2 = yAL{end}(end/2+1:end);
    % X1 = XAL{end}(1:end/2);
    % X2 = XAL{end}(end/2+1:end);
    % S = diag(Sigma);
    % S1 = S(1:end/2);
    % S2 = S(end/2+1:end);
    % % figure
    % % plot(X1,y1,'ko')
    % % hold on
    % % plot(X2,y2,'ko')
    % % plot(XNew,0*XNew,'kx','MarkerSize',12)
    % % plot(XTemp,0*XTemp,'rx','MarkerSize',12)
    % % plot(tt,f1,'b-x')
    % % plot(tt,f2,'b-+')
    % % plot(XStar1,fStar1,'r--')
    % % plot(XStar1,fStar1+2*sqrt(S1),'g--')
    % % plot(XStar1,fStar1-2*sqrt(S1),'g--')
    % % plot(XStar2,fStar2,'r--')
    % % plot(XStar2,fStar2+2*sqrt(S2),'g--')
    % % plot(XStar2,fStar2-2*sqrt(S2),'g--')

    Error(iAL) = sum((fStar(:)-ff(:)).^2);
    Error(iAL) = sqrt(Error(iAL)/(NStream*NObs));

end

if strcmp(CovType,'CovOneShare')
    fTrue = ff;
    XSample = XAL;
    ySample = yAL;
else
    fTrue = [];
    XSample = [];
    ySample = [];
end

end

%%

function [Error, tEnd] = OfflineLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound,NStream)

% CovType = [CovType 'Offline'];

% NDomain = length(X);
% rhoSource = theta0(1:2*(NDomain-1));
% theta0(1:2*(NDomain-1)) = [];
% 
% rhoTarget = theta0(1:2*(NDomain));
% theta0(1:2*(NDomain)) = [];
% 
% lambdaSource = theta0(1:2*(NDomain-1));
% theta0(1:2*(NDomain-1)) = [];
% 
% lambdaTarget = theta0(1:2*(NDomain));
% theta0(1:2*(NDomain)) = [];
% 
% sigmaSource = theta0(1:1*(NDomain-1));
% theta0(1:1*(NDomain-1)) = [];
% sigmaTarget = theta0;
% 
% % theta0 = [rhoSource(:)' rhoTarget(:)' exp(lambdaSource(:)') ...
% %      exp(lambdaTarget(:)') sigmaSource(:)' sigmaTarget(:)'];
% 
% theta0 = [rhoSource(:)' rhoTarget(:)' (lambdaSource(:)') ...
%      (lambdaTarget(:)') sigmaSource(:)' sigmaTarget(:)'];

%% Set input location for the integration

NObs = 50;
tt = linspace(Xmin,Xmax,NObs)';
XStar = repmat(tt,NStream,1);

ff = ffun(tt);

Error = zeros(1,LearningRound);
tEnd = zeros(1,LearningRound);

for i = 0:LearningRound-1
    
    NDomain = length(y);
    tStart = tic;
    yTemp = y;
    yTemp{end} = reshape(yTemp{end},[],NStream);
    yTemp{end}(end-2*(LearningRound-1-i)+1:end,:) = [];
    yTemp{end} = yTemp{end}(:);
    XTemp = X;
    XTemp{end} = reshape(XTemp{end},[],NStream);
    XTemp{end}(end-2*(LearningRound-1-i)+1:end,:) = [];
    XTemp{end} = XTemp{end}(:);
    theta = OptJOffline(theta0,yTemp,XTemp,[],CovType,NStream);
    [~, ~, Mu, Sigma] = PostInd(theta,yTemp,XTemp,[],[],[],{XStar},NDomain,CovType,NStream);
    % [Mu, Sigma] = PostOffline(theta,yTemp,XTemp,XStar,[],CovType);

    fStar = Mu;

    Error(i+1) = sum((fStar(:)-ff(:)).^2);
    Error(i+1) = sqrt(Error(i+1)/(NStream*NObs));
    
    tEnd(i+1) = toc(tStart);
end

end

