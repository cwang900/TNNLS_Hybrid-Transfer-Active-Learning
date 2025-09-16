function [Myn, Inv_Cxx_nn, Inv_KS, Det] =...
    CovCalAssembledVV(theta,y,Xf,CovType,Aim,NStream)

NDomain = length(y);

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

SeqXSource = cellfun(@length,XSource)/NStream;
SeqXSource = repmat(SeqXSource,NStream,1);
SeqXSource = mat2cell(SeqXSource,NStream,ones(1,NDomain-1));

SeqXTarget = repmat(length(XTarget{1})/NStream,NStream,NDomain);
SeqXTarget = mat2cell(SeqXTarget,NStream,ones(1,NDomain));

XS2T = [XSource XTarget];
SeqXS2T =[SeqXSource SeqXTarget{1}];
XTarget = repmat(XTarget,1,NDomain);

thetaST = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
    'UniformOutput',false);
thetaST = [thetaST thetaTarget(end)];

% \Sigma_xx^{n,n}
KTarget = CovCal(CovFun,thetaST(end),XTarget(end),XTarget(end),...
    SeqXTarget(end),SeqXTarget(end),NStream);
KTarget = KTarget{end};

% \Gamma_xx^{<n,n}
KS2T = CovCal(CovFun,thetaST,XS2T,XTarget,...
    SeqXS2T,SeqXTarget,NStream);

[mT,nT] = size(KTarget);
% if strcmp(CovType,'CovOneShare')
%     Temp = eye(2);
%     Temp = repelem(Temp,SeqXTarget{end},SeqXTarget{end});
%     KTarget = KTarget.*Temp;
% end

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
    
    % \G_xx^{n,<n}
    Prod = @(Kii,Kij) Kij'*Kii;
    Gxx_n_n = cellfun(Prod,Inv_KS,KS2T(1:NDomain-1),'UniformOutput',false);

    % \Q_xx^{n,n}
    ProdV = @(G,Kij) G*Kij;
    Qxx_nn = cellfun(ProdV,Gxx_n_n,KS2T(1:NDomain-1),'UniformOutput',false);
    Qxx_nn = cell2mat(Qxx_nn);
    Qxx_nn = reshape(Qxx_nn,mT,nT,[]);
    Qxx_nn = sum(Qxx_nn,3);
        
    % Cov(y_n|y_<n)
    Cxx_nn = KTarget-Qxx_nn;
    Cxx_nn = Cxx_nn+1e-6*eye(mT);

    % E(y_n|y_<n)
    Prod = @(Kij,y) Kij*y;
    temp = cellfun(Prod,Gxx_n_n,y(1:NDomain-1),'UniformOutput',false);
    temp = cell2mat(temp);
    Myn = sum(temp,2);
    
    % inv(Cov(y_n|y_<n))
    L_Cxx_nn = chol(Cxx_nn,'lower');
    Inv_L_Cxx_nn = L_Cxx_nn\eye(mT);
    Inv_Cxx_nn = Inv_L_Cxx_nn'*Inv_L_Cxx_nn; 
    
    if strcmp(Aim,'Optimization')
        % inv(\Sigma_xx^{<n,<n})
        KSL = cellfun(@chol,KSource,'UniformOutput',false);
        % det(\Sigma_xx^{<n,<n})
        DetFun = @(A) 2*sum(log(diag(A)));
        DetS = cellfun(DetFun,KSL);
        DetS = sum(DetS);
        % det(\Sigma_xx^{n,n})
        Det_KT = 2*sum(log(diag(L_Cxx_nn)));
        Det = DetS+Det_KT;
        Inv_Cxx_nn = Inv_L_Cxx_nn;
    end

else
    
    Inv_KS = [];
    Myn = [];
    Gxx_n_n = 0;
    Cxx_nn = KTarget;
    
    try
        L_Cxx_nn = chol(Cxx_nn,'lower');
    catch
        aa;
    end

    Inv_L_Cxx_nn = L_Cxx_nn\eye(mT);
    if strcmp(Aim,'Optimization')
        Inv_Cxx_nn = Inv_L_Cxx_nn;
        Det = 2*sum(log(diag(L_Cxx_nn)));
    else
        Inv_Cxx_nn = Inv_L_Cxx_nn'*Inv_L_Cxx_nn; %(KTT-KTS/KSS*KST)^-1
    end

end

end

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2,NStream)
    NDomain = length(theta);
    NStream = repmat({NStream},1,NDomain);
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,NStream,'UniformOutput',false);
end

    
