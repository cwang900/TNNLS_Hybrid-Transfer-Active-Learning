function [KTarget, KSTInd, KTTInd, Gxx_n_n, Inv_KS, A] =...
    CovCalInd(theta,y,Xf,XStar,CovType,NDomain,NStream)

% NDomain = length(y);

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

XSource = Xf(1:NDomain-1);
XTarget = Xf(NDomain);
XTInd = XStar;

SeqXSource = cellfun(@length,XSource)/NStream;
SeqXSource = repmat(SeqXSource,NStream,1);
SeqXSource = mat2cell(SeqXSource,NStream,ones(1,NDomain-1));

SeqXTarget = repmat(length(XTarget{1})/NStream,NStream,NDomain);
SeqXTarget = mat2cell(SeqXTarget,NStream,ones(1,NDomain));

SeqXTStar = repmat(length(XTInd{1})/NStream,NStream,NDomain);
SeqXTStar = mat2cell(SeqXTStar,NStream,ones(1,NDomain));

XS2TInd = [XSource XTInd];
SeqXS2TInd =[SeqXSource SeqXTStar{1}];

XS2Target = [XSource XTarget];
SeqXS2Target =[SeqXSource SeqXTarget{1}];

XTarget = repmat(XTarget,1,NDomain);
XTInd = repmat(XTInd,1,NDomain);

thetaST = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
    'UniformOutput',false);
thetaST = [thetaST thetaTarget(end)];

% \Sigma_xx^{n,n}
KTarget = CovCal(CovFun,thetaST(end),XTarget(end),XTarget(end),...
    SeqXTarget(end),SeqXTarget(end),NStream);
KTarget = KTarget{end};

% \Gamma_sx^{n,<n} \Gamma_ss^{n,n}
KSTInd = CovCal(CovFun,thetaST,XS2TInd,XTInd,...
    SeqXS2TInd,SeqXTStar,NStream);

% \Gamma_sx^{n,n}
KTTInd = CovCal(CovFun,thetaST(end),XTarget(end),XTInd(end),...
    SeqXTarget(end),SeqXTStar(end),NStream);
KTTInd = KTTInd{end};

% \Gamma_xx^{n,<n}
KST = CovCal(CovFun,thetaST(1:end-1),XS2Target(1:end-1),XTarget(1:end-1),...
    SeqXS2Target(1:end-1),SeqXTarget(1:end-1),NStream);

[mT,nT] = size(KTarget);
[mTS,nTS] = size(KSTInd{end});
S4Minus = sigmaTarget(1)^2*eye(nTS);
KSTInd{end} = KSTInd{end}-S4Minus;

if strcmp(CovType,'CovOneShare')
    Temp = eye(NStream);
    Temp = repelem(Temp,SeqXTarget{end},SeqXTarget{end});
    KTarget = KTarget.*Temp;
    Temp = eye(NStream);
    Temp = repelem(Temp,SeqXS2TInd{end},SeqXTStar{end});
    KSTInd{end} = KSTInd{end}.*Temp;
    Temp = eye(NStream);
    Temp = repelem(Temp,SeqXTarget{end},SeqXTStar{end});
    KTTInd = KTTInd.*Temp;
    % Temp = eye(2);
    % Temp = repelem(Temp,SeqXS2Target{end},SeqXTarget{end});
    % KST{end} = KST{end}.*Temp;
end

%%

% For transfer learning based methods
if ~isempty(thetaSource)
    % \Sigma_xx^{<n,<n}
    KSource = CovCal(CovFun,thetaSource,XSource,XSource,...
        SeqXSource,SeqXSource,NStream);
    % inv(\Sigma_xx^{<n,<n})
    Inv_KS = cellfun(@inv,KSource,'UniformOutput',false);

    % \Sigma_xx^{n,n}
    XTT = XTarget(1:end-1);
    thetaTT = thetaTarget(1:NDomain-1);
    SeqXTT = SeqXTarget(1:end-1);
    KTargetAdd = CovCal(CovFun,thetaTT,XTT,XTT,...
        SeqXTT,SeqXTT,NStream);
    KTargetAdd = cell2mat(KTargetAdd);
    KTargetAdd = reshape(KTargetAdd,mT,nT,[]);
    S4Minus = sigmaTarget(1)^2*eye(mT);
    KTargetAdd = KTargetAdd-S4Minus;
    KTargetAdd = sum(KTargetAdd,3);
    KTarget = KTarget+KTargetAdd;
    
    % \Gamma_ss^{n,n}
    XTTStar = XTInd(1:end-1);
    thetaTT = thetaTarget(1:NDomain-1);
    SeqXTTStar = SeqXTStar(1:end-1);
    KTStarAdd = CovCal(CovFun,thetaTT,XTTStar,XTTStar,...
        SeqXTTStar,SeqXTTStar,NStream);
    KTStarAdd = cell2mat(KTStarAdd);
    KTStarAdd = reshape(KTStarAdd,mTS,nTS,[]);
    S4Minus = sigmaTarget(1)^2*eye(nTS);
    KTStarAdd = KTStarAdd-S4Minus;
    KTStarAdd = sum(KTStarAdd,3);
    KSTInd{end} = KSTInd{end}+KTStarAdd;

    % \Gamma_sx^{n,n}
    KTTStarAdd = CovCal(CovFun,thetaTT,XTarget(1:end-1),XTInd(1:end-1),...
        SeqXTarget(1:end-1),SeqXTStar(1:end-1),NStream);
    KTTStarAdd = cell2mat(KTTStarAdd);
    KTTStarAdd = reshape(KTTStarAdd,mT,mTS,[]);
    KTTStarAdd = sum(KTTStarAdd,3);
    KTTInd = KTTInd+KTTStarAdd;

    % \G_xx^{n,<n}
    Prod = @(Kii,Kij) Kij'*Kii;
    Gxx_n_n = cellfun(Prod,Inv_KS,KST(1:NDomain-1),'UniformOutput',false);

    % \Q_xx^{n,n}
    ProdV = @(G,Kij) G*Kij;
    Qxx_nn = cellfun(ProdV,Gxx_n_n,KST(1:NDomain-1),'UniformOutput',false);
    Qxx_nn = cell2mat(Qxx_nn);
    Qxx_nn = reshape(Qxx_nn,mT,nT,[]);
    Qxx_nn = sum(Qxx_nn,3);
    
    % \Sigma_xx^{n,n}-Q_xx^{n,n}
    A = KTarget-Qxx_nn;

else
    
    Inv_KS = cell(1,0);
    KSTInd{end} = KSTInd{end}+1e-6*eye(mTS);
    Gxx_n_n = KTTInd/KSTInd{end};
    Gxx_n_n = {Gxx_n_n};
    A = KTarget;
    
end

end

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2,NStream)
    NDomain = length(theta);
    NStream = repmat({NStream},1,NDomain);
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,NStream,'UniformOutput',false);
end

    
