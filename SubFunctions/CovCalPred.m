function [Gxs_nn, Qxx_nn, KTStar] =...
    CovCalPred(theta,XInd,XStar,NDomain,CovType,NStream)

CovFun = str2func(CovType);

%% Prepare model parameters
if ~strcmp(CovType,'CovLMC') || length(theta) == 50
    rhoSource = reshape(theta(1:NStream*(NDomain-1)),NStream,[]);
    rhoSource = repelem(rhoSource,2,1);
    rhoSource(2:2:2*NStream,:) = 0;
    theta(1:NStream*(NDomain-1)) = [];

    rhoTarget = reshape(theta(1:NStream*(NDomain)),NStream,[]);
    rhoTarget = repelem(rhoTarget,2,1);
    rhoTarget(2:2:2*NStream,:) = 0;
    theta(1:NStream*(NDomain)) = [];

    lambdaSource = reshape(theta(1:NStream*(NDomain-1)),NStream,[]);
    lambdaSource = repelem(lambdaSource,2,1);
    theta(1:NStream*(NDomain-1)) = [];

    lambdaTarget = reshape(theta(1:NStream*(NDomain)),NStream,[]);
    lambdaTarget = repelem(lambdaTarget,2,1);
    theta(1:NStream*(NDomain)) = [];

    sigmaSource = reshape(theta(1:1*(NDomain-1)),1,[]);
    sigmaSource = repmat(sigmaSource,NStream,1);
    % if NDomain == 1
    %     sigmaSource = repmat(theta(1:1*(NDomain-1)),2,1);
    % else
    %     sigmaSource = repmat(theta(1:1*(NDomain-1))',2,1);
    % end
    theta(1:1*(NDomain-1)) = [];
    sigmaTarget = repmat(theta,NStream,NDomain);
elseif length(theta) == 58
    rhoSource = reshape(theta(1:NStream*(NDomain-1)),NStream,[]);
    rhoSource = repelem(rhoSource,2,1);
    rhoSource(2:2:2*NStream,:) = 0;
    theta(1:NStream*(NDomain-1)) = [];

    rhoTarget = reshape(theta(1:NStream*(NDomain+1)),NStream,[]);
    rhoTarget = repelem(rhoTarget,2,1);
    rhoTarget(2:2:2*NStream,:) = 0;
    rhoTarget(2:2:2*NStream,end-1) = rhoTarget(1:2:2*NStream,end);
    rhoTarget(:,end) = [];
    theta(1:NStream*(NDomain+1)) = [];

    lambdaSource = reshape(theta(1:NStream*(NDomain-1)),NStream,[]);
    lambdaSource = repelem(lambdaSource,2,1);
    theta(1:NStream*(NDomain-1)) = [];

    lambdaTarget = reshape(theta(1:NStream*(NDomain)),NStream,[]);
    lambdaTarget = repmat(lambdaTarget,2,1);
    theta(1:NStream*(NDomain)) = [];

    sigmaSource = reshape(theta(1:1*(NDomain-1)),1,[]);
    sigmaSource = repmat(sigmaSource,NStream,1);

    theta(1:1*(NDomain-1)) = [];
    sigmaTarget = repmat(theta,NStream,NDomain);
else
    rhoSource = reshape(theta(1:NStream*(NDomain-1)),NStream,[]);
    rhoSource = repelem(rhoSource,2,1);
    rhoSource(2:2:2*NStream,:) = 0;
    theta(1:NStream*(NDomain-1)) = [];

    % rhoTarget = reshape(theta(1:NStream*(NDomain+1)),NStream,[]);
    % rhoTarget = repelem(rhoTarget,2,1);
    % rhoTarget(2:2:2*NStream,:) = 0;
    % rhoTarget(2:2:2*NStream,end-1) = rhoTarget(1:2:2*NStream,end);
    % rhoTarget(:,end) = [];
    % theta(1:NStream*(NDomain+1)) = [];

    rhoTarget = reshape(theta(1:NStream*(NDomain)),NStream,[]);
    rhoTarget = repelem(rhoTarget,2,1);
    rhoTarget(2:2:2*NStream,:) = 0;
    theta(1:NStream*(NDomain)) = [];
    
    lambdaSource = reshape(theta(1:NStream*(NDomain-1)),NStream,[]);
    lambdaSource = repelem(lambdaSource,2,1);
    theta(1:NStream*(NDomain-1)) = [];

    lambdaTarget = reshape(theta(1:NStream*(NDomain)),NStream,[]);
    lambdaTarget = repmat(lambdaTarget,2,1);
    theta(1:NStream*(NDomain)) = [];

    sigmaSource = reshape(theta(1:1*(NDomain-1)),1,[]);
    sigmaSource = repmat(sigmaSource,NStream,1);

    theta(1:1*(NDomain-1)) = [];
    sigmaTarget = repmat(theta,NStream,NDomain);
end

thetaSource = [rhoSource;lambdaSource;sigmaSource];
thetaSource = mat2cell(thetaSource,5*NStream,ones(1,NDomain-1));
thetaTarget = [rhoTarget;lambdaTarget;sigmaTarget];
thetaTarget = mat2cell(thetaTarget,5*NStream,ones(1,NDomain));

%%

XTInd = XInd;
XTStar = XStar;

SeqXInd = repmat(length(XTInd{1})/NStream,NStream,NDomain);
SeqXInd = mat2cell(SeqXInd,NStream,ones(1,NDomain));

SeqXTStar = repmat(length(XTStar{1})/NStream,NStream,NDomain);
SeqXTStar = mat2cell(SeqXTStar,NStream,ones(1,NDomain));

XTInd = repmat(XTInd,1,NDomain);
XTStar = repmat(XTStar,1,NDomain);

% \Gamma_ss_nn
KTInd = CovCal(CovFun,thetaTarget(end),XTInd(end),XTInd(end),...
    SeqXInd(end),SeqXInd(end),NStream);
KTInd = KTInd{end};

% \Gamma_sx_nn
KTIndTStar = CovCal(CovFun,thetaTarget(end),XTInd(end),XTStar(end),...
    SeqXInd(end),SeqXTStar(end),NStream);
KTIndTStar = KTIndTStar{end};

% \Gamma_xx_nn
KTStar = CovCal(CovFun,thetaTarget(end),XTStar(end),XTStar(end),...
    SeqXTStar(end),SeqXTStar(end),NStream);
KTStar = KTStar{end};

[mT,nT] = size(KTInd);
[mTS,nTS] = size(KTStar);
S4Minus = sigmaTarget(1)^2*eye(mT);
KTInd = KTInd-S4Minus;
if strcmp(CovType,'CovOneShare')
    Temp = eye(NStream);
    Temp = repelem(Temp,SeqXInd{end},SeqXInd{end});
    KTInd = KTInd.*Temp;
    Temp = eye(NStream);
    Temp = repelem(Temp,SeqXTStar{end},SeqXTStar{end});
    KTStar = KTStar.*Temp;
    Temp = eye(NStream);
    Temp = repelem(Temp,SeqXInd{end},SeqXTStar{end});
    KTIndTStar = KTIndTStar.*Temp;
end

%%

% For transfer learning based methods
if ~isempty(thetaSource)

    % \Gamma_ss_nn
    XTT = XTInd(1:end-1);
    thetaTT = thetaTarget(1:NDomain-1);
    SeqXTT = SeqXInd(1:end-1);
    KTIndAdd = CovCal(CovFun,thetaTT,XTT,XTT,...
        SeqXTT,SeqXTT,NStream);
    KTIndAdd = cell2mat(KTIndAdd);
    KTIndAdd = reshape(KTIndAdd,mT,nT,[]);
    S4Minus = sigmaTarget(1)^2*eye(mT);
    KTIndAdd = KTIndAdd-S4Minus;
    KTIndAdd = sum(KTIndAdd,3);
    KTInd = KTInd+KTIndAdd;

    % \Gamma_xx_nn
    XTTStar = XTStar(1:end-1);
    thetaTT = thetaTarget(1:NDomain-1);
    SeqXTTStar = SeqXTStar(1:end-1);
    KTStarAdd = CovCal(CovFun,thetaTT,XTTStar,XTTStar,...
        SeqXTTStar,SeqXTTStar,NStream);
    KTStarAdd = cell2mat(KTStarAdd);
    KTStarAdd = reshape(KTStarAdd,mTS,nTS,[]);
    S4Minus = sigmaTarget(1)^2*eye(nTS);
    KTStarAdd = KTStarAdd-S4Minus;
    KTStarAdd = sum(KTStarAdd,3);
    KTStar = KTStar+KTStarAdd-S4Minus;

    % \Gamma_sx_nn
    KTIndTStarAdd = CovCal(CovFun,thetaTT,XTInd(1:end-1),XTStar(1:end-1),...
        SeqXInd(1:end-1),SeqXTStar(1:end-1),NStream);
    KTIndTStarAdd = cell2mat(KTIndTStarAdd);
    KTIndTStarAdd = reshape(KTIndTStarAdd,mT,mTS,[]);
    KTIndTStarAdd = sum(KTIndTStarAdd,3);
    KTIndTStar = KTIndTStar+KTIndTStarAdd;
    
end
KTInd = KTInd+1e-6*eye(mT);
Gxs_nn = KTIndTStar'/KTInd;
Qxx_nn = Gxs_nn*KTIndTStar;

end

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2,NStream)
    NDomain = length(theta);
    NStream = repmat({NStream},1,NDomain);
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,NStream,'UniformOutput',false);
end


    
