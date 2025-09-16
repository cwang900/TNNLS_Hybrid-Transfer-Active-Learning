clear
clc
% close all

%%

NRep = 100;
LearningRound = 26;
NoiseSigma = 1;
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
% load theta0

parfor iRep = 1:NRep
    %% Generate data

    NDomain = 8;
    NObs = 30;
    Xmin = 0;
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

    g1 = @(b,x) b(1)+b(2)*x+b(3)*x.^2+b(4)*x.^3;
    g2 = @(w,rho,x) 2*(rho*sin(w*x)+(1-rho)*cos(w*x+pi/3));
    g = @(rho,x) 0*x.^3+x.^2-x+2+4*(rho*sin(1*x)+(1-rho)*cos(1*x));

    Fun = {g1,g2};

    % b = [5 -5 4 -3 0 0 0 2; 2 1 0 0 2 3 2 -1; 1 -2 -1 1/2 0 0 1 -1;...
    %      0 0 3 -1 -3 2 2 1];%
    b = [unifrnd(0,5,1,NDomain); unifrnd(0,2,3,NDomain);];%
    b(4,:) = 0;
    Ind = randi([1 4],1,NDomain);
    Ind = sub2ind(size(b),Ind,1:NDomain);
    b(Ind) = 0;
    % w = normrnd(2.5,0.2,1,NDomain);
    w = normrnd(1,0.1,1,NDomain);
    rho = [1 0; 0.2 0.8];%
    
    f = cell(1,NDomain);
    % figure(1)
    % hold on
    % figure(2)
    % hold on
    for i = 1:NDomain
        X{i} = [t{i};t{i}];
        if i<NDomain
            f{i} = [g1(b(:,i),t{i})+g2(w(i),rho(1,1),t{i});...
                g1([-b(1,i); b(2:end,i)],t{i})+g2(w(i),rho(1,2),t{i})];
        else
            f{i} = [g(rho(2,1),t{i}); g(rho(2,2),t{i})];
        end
        y{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
        % figure(1)
        % plot(t{i},f{i}(1:end/2))
        % figure(2)
        % plot(t{i},f{i}(end/2+1:end))
    end

    ffun = @(t) [g(rho(2,1),t) ...
        g(rho(2,2),t)];

    %%

    XRand = unifrnd(Xmin,Xmax,2,LearningRound);
   
    CovType = {'CovOneShare','CovSGPTransfer','CovNonTransfer','CovLMC'};
    % FunVar{5} = CovType;

    %% Set intial parameter

    NDomain = length(X);
    Rho0 = [unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = 1*([unifrnd(5,5,1,NDomain-1);unifrnd(5,5,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(2,2,2,NDomain-1) unifrnd(2,2,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%thetaSample;% 

    %% Active learning

    try
    [Error(iRep,:), fTrue, XSample, ySample, tEndProp(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{1},LearningRound,XRand,NoiseSigma);

    NDomain = length(X);
    Rho0 = [unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(5,5,1,NDomain-1);unifrnd(4,4,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)]/sqrt(NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = 1*([unifrnd(5,5,2,NDomain-1) unifrnd(2,2,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 0.1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];

    [ErrorSGPTransfer(iRep,:), ~, ~, ~, tEndSGPTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{2},LearningRound,XRand,NoiseSigma);

    NDomain = length(X);
    Rho0 = [unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(3,3,1,NDomain-1);unifrnd(4,4,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)]/sqrt(NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([0.5*unifrnd(3,3,2,NDomain-1) 0.2*unifrnd(4,4,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 0.1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(1,:)'; Scale(end); sigma(:)];
    [ErrorOffline(iRep,:), tEnd(iRep,:)] = OfflineLearning...
        (XSample,ySample,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound);

    [ErrorLMC(iRep,:), ~, ~, ~, tEndLMC(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,XRand,NoiseSigma);

    NDomain = length(X);
    Rho0 = [unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(3,3,1,NDomain-1);unifrnd(4,4,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)]/sqrt(NDomain);%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([0.5*unifrnd(3,3,2,NDomain-1) 0.5*unifrnd(4,4,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 0.1*ones(1,NDomain);%

    y = y(end);
    X = X(end);
    Rho = Rho(:,1);%RhoSample(:,end);%unifrnd(1,2,2,NDomain);%
    Scale = Scale(:,end);%unifrnd(0.1,0.2,2,1);%ScaleSample(:,end);%
    sigma = sigma(1);
    theta0 = [Rho(:); Scale(:); sigma(:)];
    [ErrorNonTransfer(iRep,:), ~, ~, ~, tEndNonTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{3},LearningRound,XRand,NoiseSigma);

    catch
        continue;
    end
end
save('RandomPoly')

%%
function [Error, fTrue, XSample, ySample, tEnd] = ActiveLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound,XRand,NoiseSigma)
tEnd = zeros(1,LearningRound);
NDomain = length(X);
% b = b(:,NDomain:NDomain:end);

%% Set input location for the integration

NObs = 50;
tt = linspace(Xmin,Xmax,NObs)';
XStar = tt;
% XStar = sort(unifrnd(Xmax/2,Xmax,NObs/2-2,1));
XStar = [XStar;XStar];

NInd = 40;
XTInd = linspace(Xmin,Xmax,NInd)';
XTInd = [XTInd;XTInd];

%% Parameter optimization & Prediction

theta = OptJOffline(theta0,y,X,XTInd,CovType);

alpha = [];
C = [];
alphat = zeros(2*NInd,LearningRound);
Ct = zeros(2*NInd,2*NInd,LearningRound);
[Mu, Sigma, alphat(:,1), Ct(:,:,1)] = PostInd(theta,y,X,{XStar},alpha,C,{XTInd},NDomain,CovType);

%%

fStar1 = Mu(1:end/2);
fStar2 = Mu(end/2+1:end);

ff = ffun(tt);
f1 = ff(:,1);
f2 = ff(:,2);
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
% plot(tt,f1,'b-x')
% plot(tt,f2,'b-+')
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
Error = zeros(1,LearningRound);
Error(1) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
Error(1) = sqrt(Error(1)/(2*NObs));

while iAL<LearningRound %Error > 1e-3
    tstart = tic;
    p = sort(lhsdesign(NSample,2),2);
    XSample0 = Xmin+p*(Xmax-Xmin);
    % XTemp = OPTIMSE(XSample0,XStar,Ct(:,:,iAL),XTInd,theta,NDomain,CovType,Xmin,Xmax);
    XTemp = XRand(:,iAL);
    XNew = [XNew XTemp(:)];
    ffNew = ffun(XTemp(:));
    yNew = ffNew;%+normrnd(0,thetaAL(end,iAL),1,2);
    yNew = yNew+normrnd(0,NoiseSigma,size(yNew));
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
            PostInd(theta,yAL,XAL,{XStar},alphat(:,iAL),Ct(:,:,iAL),{XTInd},NDomain,CovType);
    else
        [Mu, Sigma, alphat(:,iAL+1), Ct(:,:,iAL+1)] = ...
            PostInd(theta,{yt},{Xt},{XStar},alphat(:,iAL),Ct(:,:,iAL),{XTInd},NDomain,CovType);
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
    % plot(tt,f1,'b-x')
    % plot(tt,f2,'b-+')
    % plot(XStar1,fStar1,'r--')
    % plot(XStar1,fStar1+2*sqrt(S1),'g--')
    % plot(XStar1,fStar1-2*sqrt(S1),'g--')
    % plot(XStar2,fStar2,'r--')
    % plot(XStar2,fStar2+2*sqrt(S2),'g--')
    % plot(XStar2,fStar2-2*sqrt(S2),'g--')

    Error(iAL) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
    Error(iAL) = sqrt(Error(iAL)/(2*NObs));

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

NObs = 50;
tt = linspace(Xmin,Xmax,NObs)';
XStar = tt;
% XStar = sort(unifrnd(Xmax/2,Xmax,NObs/2-2,1));
XStar = [XStar;XStar];

ff = ffun(tt);
f1 = ff(:,1);
f2 = ff(:,2);

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
