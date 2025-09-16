clear
clc

%%

NRep = 100;
LearningRound = 30;
NoiseSigma = 0.5;
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
FunVar = cell(1,5);

% parpool(10) % Set parallel computing
parfor iRep = 1:NRep
    %% Generate data

    NDomain = 3;
    NStream = 2;
    NObs = 30;
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

    % g1 = @(u,x,rho) u+2*cos(x);%b(1)+b(2)*x.^1;
    % g2 = @(u,x,rho) u+3*cos(x);
    % g = @(u,x,eta,rho) u+(2*cos(x));

    % g1 = @(u,x,rho) u+2*sin(x)+rho*sin(2*x)+(1-rho)*sin(2*x+pi/2);%b(1)+b(2)*x.^1;
    % g2 = @(u,x,rho) u+2*cos(x)+rho*sin(3*x)+(1-rho)*sin(3*x+pi/3);
    % g = @(u,x,eta,rho) u+exp(-x./eta).*(rho*5*sin(x)+(1-rho)*5*cos(x));

    % g1 = @(u,x,rho) u+2*sin(x)+rho*sin(2*x)+(1-rho)*sin(2*x+pi/2);
    % g2 = @(u,x,rho) u+2*cos(x)+rho*sin(3*x)+(1-rho)*sin(3*x+pi/3);
    % g = @(u,x,Phi,eta,rho) u+3*sin(2.5*x)+exp(-x./eta).*(rho*5*sin(x+Phi)+(1-rho)*5*cos(x));
    
    w0 = 0.7;
    g1 = @(u,x,rho) u+2*((1-rho)*sin(w0*x)+rho*sin(1*x));
    g2 = @(u,x,rho) u+2*((1-rho)*cos(w0*x)+rho*cos(1*x));
    g = @(u,x,Phi,eta,rho) u+...
        exp(-x./eta).*(rho*2*sin(w0*x)+(1-rho)*2*cos(w0*x));

    Fun = {g1,g2};

    b = {[2 1; 0 1]; [1.8 1; -0.2 1]; [1.6 1; -0.1 0.9]};%[ones(NDomain,1) unifrnd(1/3,1/2,NDomain,1)];
    w = [3 2.5; 2.6 2; 2.4 2.2];%normrnd(2.5,0.3,NDomain,2);%[0.9 1.3; 1.1 1.2; 1 1.1];%unifrnd(1,2,NDomain,2);1.2*ones(NDomain,1);%
    Phi = [0,3*pi/4];
    % rho = [[0.5; 0.7] [0.6; 0.5] [0.4; 0.8]];%[0.8 0.2 0.6];%
    % eta = [5,-8];
    rho = [[0.5; 0.7] [0.5; 0.7] [0.7; 0.3]];
    eta = [5,-10];
    U = [unifrnd(0,5,1,NDomain);unifrnd(-5,0,1,NDomain);];%zeros(2,NDomain);%
    
    f = cell(1,NDomain);
    for i = 1:NDomain
        X{i} = [t{i};t{i}];
        if i<NDomain
            f{i} = [Fun{i}(U(1,i),t{i},rho(1,i)); Fun{i}(U(2,i),t{i},rho(2,i))];
        else
            f{i} = [g(U(1,i),t{i},Phi(1),eta(1),rho(1,i)); g(U(2,i),t{i},Phi(2),eta(2),rho(2,i))];
        end
        y{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
    end

    ffun = @(t) [g(U(1,NDomain),t,Phi(1),eta(1),rho(1,NDomain)) ...
        g(U(2,NDomain),t,Phi(2),eta(2),rho(2,NDomain))];

    %%

    XRand = unifrnd(Xmin,Xmax,2,LearningRound);
   
    CovType = {'CovOneShare','CovSGPTransfer','CovNonTransfer','CovLMC'};
    % FunVar{5} = CovType;

    %% Set intial parameter
    
    NDomain = length(X);
    Rho0 = 1*[unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = [1*unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(5,5,2,NDomain-1) unifrnd(2,2,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 1*ones(1,NDomain);%

    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];%thetaSample;% 

    %% Active learning
    
    try
    [RMSE(iRep,:), MAPE(iRep,:), fTrue, XSample, ySample] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{1},LearningRound,XRand,NoiseSigma,NStream);

    NDomain = length(X);
    Rho0 = 1*[unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*[unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(5,5,2,NDomain-1) unifrnd(5,5,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0(:); Scale(:); sigma(:)];

    [RMSESGPTransfer(iRep,:), MAPESGPTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{2},LearningRound,XRand,NoiseSigma,NStream);
    
    NDomain = length(X);
    Rho0 = 1*[unifrnd(1,1.5,1,NDomain-1);unifrnd(1,1.5,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*[unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(0.5,0.5,2,NDomain-1) unifrnd(0.5,0.5,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 1*ones(1,NDomain);%
    theta0 = [Rho0(:); Rho(:); Scale0(1,:)'; Scale(end); sigma(:)];
    [RMSELMC(iRep,:), MAPELMC(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{4},LearningRound,XRand,NoiseSigma,NStream);

    NDomain = length(X);
    Rho0 = 1*[unifrnd(1,1,1,NDomain-1);unifrnd(1,1,1,NDomain-1)];%unifrnd(1,3,2,NDomain-1);%
    Scale0 = ([unifrnd(4,4,1,NDomain-1);unifrnd(4,4,1,NDomain-1)]);%unifrnd(2,3,2,NDomain-1);%
    Rho = 1*[unifrnd(1,1,2,NDomain-1) 1*unifrnd(1,1,2,1)];%[lmin+dl*Corr unifrnd(1,2,2,1)]/((NDomain).^(1/3));%unifr  nd(1,2,2,NDomain);%
    Scale = ([unifrnd(1,1,2,NDomain-1) unifrnd(1,1,2,1);]);%log([2-Corr unifrnd(1,2,2,1)]);%unifrnd(5,6,1,NDomain)]);
    beta = 0.5;
    sigma = 1*ones(1,NDomain);%

    y = y(end);
    X = X(end);
    Rho = Rho(:,end);%RhoSample(:,end);%unifrnd(1,2,2,NDomain);%
    Scale = Scale(:,end);%unifrnd(0.1,0.2,2,1);%ScaleSample(:,end);%
    sigma = sigma(1);
    theta0 = [Rho(:); Scale(:); sigma(:)];
    [RMSENonTransfer(iRep,:), MAPENonTransfer(iRep,:)] = ActiveLearning...
        (X,y,ffun,Xmin,Xmax,theta0,CovType{3},LearningRound,XRand,NoiseSigma,NStream);

    catch
        continue;
    end
end
save('ALTrig')

%%
function [RMSE, MAPE, fTrue, XSample, ySample] = ActiveLearning(X,y,ffun,Xmin,Xmax,theta0,CovType,LearningRound,XRand,NoiseSigma,NStream)

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

theta = OptJOffline(theta0,y,X,XTInd,CovType,NStream);

alpha = [];
C = [];
alphat = zeros(2*NInd,LearningRound);
Ct = zeros(2*NInd,2*NInd,LearningRound);
[Mu, Sigma, alphat(:,1), Ct(:,:,1)] = PostInd(theta,y,X,{XStar},alpha,C,{XTInd},NDomain,CovType,NStream);

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
RMSE = zeros(1,LearningRound);
RMSE(1) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
RMSE(1) = sqrt(RMSE(1)/(2*NObs));

MAPE = zeros(1,LearningRound);
MAPE(1) = sum(abs((fStar1(:)-f1(:))./f1(:)))+sum(abs((fStar2(:)-f2(:))./f2(:)));
MAPE(1) = 100*(MAPE(1)/(2*NObs));

while iAL<LearningRound %Error > 1e-3

    p = sort(lhsdesign(NSample,2),2);
    % p = sobolset(2,'Skip',1e2,'Leap',1e2);
    % p = sort(net(p,1),2);
    XSample0 = Xmin+p*(Xmax-Xmin);
    % XTemp = OPTIMSE(XSample0,XStar,Ct(:,:,iAL),XTInd,theta,NDomain,CovType,Xmin,Xmax,NStream);
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
    % plot(XTemp,0*XTemp,'rx','MarkerSize',12)
    % plot(tt,f1,'b-x')
    % plot(tt,f2,'b-+')
    % plot(XStar1,fStar1,'r--')
    % plot(XStar1,fStar1+2*sqrt(S1),'g--')
    % plot(XStar1,fStar1-2*sqrt(S1),'g--')
    % plot(XStar2,fStar2,'r--')
    % plot(XStar2,fStar2+2*sqrt(S2),'g--')
    % plot(XStar2,fStar2-2*sqrt(S2),'g--')

    RMSE(iAL) = sum((fStar1(:)-f1(:)).^2)+sum((fStar2(:)-f2(:)).^2);
    RMSE(iAL) = sqrt(RMSE(iAL)/(2*NObs));

    MAPE(iAL) = sum(abs((fStar1(:)-f1(:))./f1(:)))+sum(abs((fStar2(:)-f2(:))./f2(:)));
    MAPE(iAL) = 100*(MAPE(iAL)/(2*NObs));

end

if strcmp(CovType,'CovOneShare')
    fTrue = [f1 f2];
    XSample = XAL;
    ySample = yAL;
end

end
