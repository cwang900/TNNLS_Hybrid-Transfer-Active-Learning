function [Mu, Sigma] = PostOffline2D(theta,y,Xf,XStar,Xu,CovType)

NDomain = length(y);
if strcmp(CovType,'CovLMC')
    theta = theta';
    rhoSource = theta(1:2*(NDomain-1));
    theta(1:2*(NDomain-1)) = [];

    rhoTarget = theta(1:2*(NDomain));
    theta(1:2*(NDomain)) = [];

    lambdaSource = theta(1:NDomain-1);
    lambdaSource = repelem(lambdaSource,2,1);
    theta(1:NDomain-1) = [];

    lambdaTarget = theta(1);
    lambdaTarget = repelem(lambdaTarget,2,1);
    lambdaTarget = [lambdaSource lambdaTarget];
    theta(1) = [];

    sigma = theta;
    
    theta = [rhoSource rhoTarget lambdaSource(:)' lambdaTarget(:)' sigma];
    theta = theta(:);
end
Aim = '';
if NDomain == 1
    [~, Inv_KT_KST, ~, KTStar, KTTStar] = ...
        CovCalAssembled(theta,y,Xf,XStar,CovType,Aim);
    yKT = Inv_KT_KST*y{end};
    KTStarT = KTStar{1};
    Mu = KTStarT'*yKT;
    KKTTStar = KTStarT'*Inv_KT_KST*KTStarT;
    Sigma = KTTStar-KKTTStar;
else
    [Inv_KS_Full, Inv_KT_KST, KST_InvOT, KTStar, KTTStar, KTStarT] = ...
        CovCalAssembled(theta,y,Xf,XStar,CovType,Aim);
    yS = cell2mat(y(1:NDomain-1)');
    yKS = Inv_KS_Full*yS;
    yS_KST = @(yS,Inv_KST) Inv_KST'*yS;
    yKST = cellfun(yS_KST,...
        y(1:NDomain-1),KST_InvOT,'UniformOutput',false);
    yKST = cell2mat(yKST);
    yKST = sum(yKST,2);
    yT_KST = @(Inv_KST) Inv_KST*y{NDomain};
    yTKST = cellfun(yT_KST,...
        KST_InvOT,'UniformOutput',false);
    yTKST = cell2mat(yTKST');
    yKT = Inv_KT_KST*y{end};

    KTStarS = cell2mat(KTStar(1:NDomain-1)')';
    MuS = KTStarS*yKS;
    MuST = KTStarS*yTKST+KTStarT'*yKST;
    MuTT = KTStarT'*yKT;
    Mu = MuS+MuST+MuTT;

    KKTStarS = KTStarS*Inv_KS_Full*KTStarS';
    KKTStarT = KTStarT'*cell2mat(KST_InvOT')'*KTStarS';
    KKTStarT = KKTStarT+KKTStarT';
    KKTTStar = KTStarT'*Inv_KT_KST*KTStarT;
    Sigma = KTTStar-(KKTStarS+KKTStarT+KKTTStar);
end
end
% %% Prepare model parameters
% 
% rhoSource = reshape(theta(1:2*(NDomain-1)),2,[]);
% rhoSource = repelem(rhoSource,2,1);
% rhoSource([2 4],:) = 0;
% theta(1:2*(NDomain-1)) = [];
% 
% rhoTarget = reshape(theta(1:2*(NDomain)),2,[]);
% rhoTarget = repelem(rhoTarget,2,1);
% rhoTarget([2 4],:) = 0;
% theta(1:2*(NDomain)) = [];
% 
% lambdaSource = reshape(theta(1:2*(NDomain-1)),2,[]);
% lambdaSource = repelem(lambdaSource,2,1);
% theta(1:2*(NDomain-1)) = [];
% 
% lambdaTarget = reshape(theta(1:2*(NDomain)),2,[]);
% lambdaTarget = repelem(lambdaTarget,2,1);
% theta(1:2*(NDomain)) = [];
% 
% if NDomain == 1
%     sigmaSource = repmat(theta(1:1*(NDomain-1)),2,1);
% else
%     sigmaSource = repmat(theta(1:1*(NDomain-1))',2,1);
% end
% theta(1:1*(NDomain-1)) = [];
% sigmaTarget = repmat(theta,2,NDomain);
% 
% thetaSource = [rhoSource;lambdaSource;sigmaSource];
% thetaSource = mat2cell(thetaSource,10,ones(1,NDomain-1));
% thetaTarget = [rhoTarget;lambdaTarget;sigmaTarget];
% thetaTarget = mat2cell(thetaTarget,10,ones(1,NDomain));
% 
% %% Covariance calculation
% 
% CovFun = str2func(CovType);
% 
% XSource = Xf(1:NDomain-1);
% XTarget = Xf(NDomain);
% 
% SeqXSource = cellfun(@length,XSource)/2;
% SeqXSource = repmat(SeqXSource,2,1);
% SeqXSource = mat2cell(SeqXSource,2,ones(1,NDomain-1));
% 
% SeqXTarget = repmat(length(XTarget{1})/2,2,NDomain);
% SeqXTarget = mat2cell(SeqXTarget,2,ones(1,NDomain));
% 
% XS2T = [XSource XTarget];
% SeqXS2T =[SeqXSource SeqXTarget{1}];
% XTarget = repmat(XTarget,1,NDomain);
% 
% thetaST = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
%     'UniformOutput',false);
% thetaST = [thetaST thetaTarget(end)];
% 
% % KTarget = cellfun(CovFun,thetaST,XS2T,XTarget,...
% %     SeqXS2T,SeqXTarget,'UniformOutput',false);
% KTarget = CovCal(CovFun,thetaST,XS2T,XTarget,...
%     SeqXS2T,SeqXTarget);
% [m,n] = size(KTarget{end});
% 
% XTStar = repmat({XStar},1,NDomain);
% thetaTStar = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
%     'UniformOutput',false);
% thetaTStar = [thetaTStar, thetaTarget(end)];
% SeqXTStar = repmat(NStar/2,2,NDomain);
% SeqXTStar = mat2cell(SeqXTStar,2,ones(1,NDomain));
% 
% % KTStar = cellfun(CovFun,thetaTStar,Xf,XTStar,...
% %     SeqXS2T,SeqXTStar,'UniformOutput',false);
% KTStar = CovCal(CovFun,thetaTStar,Xf,XTStar,...
%     SeqXS2T,SeqXTStar);
% 
% % KTTStar = cellfun(CovFun,thetaTarget,XTStar,XTStar,...
% %     SeqXTStar,SeqXTStar,'UniformOutput',false);
% KTTStar = CovCal(CovFun,thetaTarget,XTStar,XTStar,...
%     SeqXTStar,SeqXTStar);
% KTTStar = cell2mat(KTTStar);
% KTTStar = reshape(KTTStar,NStar,NStar,[]);
% S4Minus = sigmaTarget(1)^2*eye(NStar);
% KTTStar = KTTStar-S4Minus;
% KTTStar = sum(KTTStar,3);
% 
% % For transfer learning based methods
% if ~isempty(thetaSource)
% 
%     % KSource = cellfun(CovFun,thetaSource,XSource,XSource,...
%     %     SeqXSource,SeqXSource,'UniformOutput',false);
%     KSource = CovCal(CovFun,thetaSource,XSource,XSource,...
%         SeqXSource,SeqXSource);    
% 
%     XTT = XTarget(1:end-1);
%     thetaTT = thetaTarget(1:NDomain-1);
%     SeqXTT = SeqXTarget(1:end-1);
% 
%     % KTargetAdditional = cellfun(CovFun,thetaTT,XTT,XTT,...
%     %     SeqXTT,SeqXTT,'UniformOutput',false);
%     KTargetAdditional = CovCal(CovFun,thetaTT,XTT,XTT,...
%         SeqXTT,SeqXTT);
%     KTargetAdditional = cell2mat(KTargetAdditional);
%     KTargetAdditional = reshape(KTargetAdditional,m,n,[]);
%     S4Minus = sigmaTarget(1)^2*eye(m);
%     KTargetAdditional = KTargetAdditional-S4Minus;
%     KTargetAdditional = sum(KTargetAdditional,3);
%     KTarget{end} = KTarget{end}+KTargetAdditional;
% 
%     % KTStarAdditional = cellfun(CovFun,thetaTarget(1:NDomain-1),...
%     %     XTT,XTStar(1:NDomain-1),...
%     %     SeqXTT,SeqXTStar(1:NDomain-1),'UniformOutput',false);
%     KTStarAdditional = CovCal(CovFun,thetaTarget(1:NDomain-1),...
%         XTT,XTStar(1:NDomain-1),SeqXTT,SeqXTStar(1:NDomain-1));
%     KTStarAdditional = cell2mat(KTStarAdditional);
%     KTStarAdditional = reshape(KTStarAdditional,sum(SeqXTT{1}),sum(SeqXTStar{1}),[]);
%     KTStarAdditional = sum(KTStarAdditional,3);
%     KTStar{end} = KTStar{end}+KTStarAdditional;
%     KTStarT = KTStar{end};
% 
%     G = cellfun(@Prod,KSource,KTarget(1:NDomain-1),'UniformOutput',false);
%     Gtemp = cell2mat(G);
%     % KSKTKS = Gtemp*cell2mat(KTarget(1:NDomain-1));
%     KST = cellfun(@ProdV,G,KTarget(1:NDomain-1),'UniformOutput',false);
%     KST = cell2mat(KST);
%     KST = reshape(KST,m,n,[]);
%     KST = sum(KST,3);
% 
%     KT_KST = KTarget{end}-KST;
%     KT_KST = KT_KST+1e-6*eye(m);
%     L_KT_KST = chol(KT_KST,'lower');
% 
%     Inv_L_KT_KST = L_KT_KST\eye(m);
%     Inv_KT_KST = Inv_L_KT_KST'*Inv_L_KT_KST; %(KTT-KTS/KSS*KST)^-1Inv_KS = cellfun(@inv,KSource,'UniformOutput',false);
% 
%     Inv_KS = cellfun(@inv,KSource,'UniformOutput',false);
%     KST_InvOT = cellfun(@(G)ProdVV(G,Inv_KT_KST),G,...
%         'UniformOutput',false);
%     KST_InvOT_mat = cell2mat(KST_InvOT');
%     Inv_KS_Sub = -KST_InvOT_mat*Gtemp;
%     Inv_KS_Full = blkdiag(Inv_KS{:})+Inv_KS_Sub;
%     yS = cell2mat(y(1:NDomain-1)');
%     yKS = Inv_KS_Full*yS;
%     yKST = cellfun(@(yS,Inv_KST)yS_KST(yS,Inv_KST),...
%         y(1:NDomain-1),KST_InvOT,'UniformOutput',false);
%     yKST = cell2mat(yKST);
%     yKST = sum(yKST,2);
%     yTKST = cellfun(@(Inv_KST)yT_KST(y{NDomain},Inv_KST),...
%         KST_InvOT,'UniformOutput',false);
%     yTKST = cell2mat(yTKST');
%     yKT = Inv_KT_KST*y{end};
% 
%     KTStarS = cell2mat(KTStar(1:NDomain-1)')';
%     MuS = KTStarS*yKS;
%     MuST = KTStarS*yTKST+KTStarT'*yKST;
%     MuTT = KTStarT'*yKT;
%     Mu = MuS+MuST+MuTT;
% 
%     KKTStarS = KTStarS*Inv_KS_Full*KTStarS';
%     KKTStarT = KTStarT'*cell2mat(KST_InvOT')'*KTStarS';
%     KKTStarT = KKTStarT+KKTStarT';
%     KKTTStar = KTStarT'*Inv_KT_KST*KTStarT;
%     Sigma = KTTStar-(KKTStarS+KKTStarT+KKTTStar);
% 
% else
% 
%     % KTarget = CovCal(CovFun,thetaTarget,XTarget,XTarget,...
%     %     SeqXTarget,SeqXTarget);
%     % 
%     % thetaTStar = thetaTarget;
%     % SeqXTStar = repmat(NStar/2,2,NDomain);
%     % SeqXTStar = mat2cell(SeqXTStar,2,ones(1,NDomain));
%     % XTStar = repmat({XStar},1,NDomain);
%     % SeqXS2T = SeqXTarget;
%     % 
%     % KTStar = CovCal(CovFun,thetaTStar,Xf,XTStar,...
%     %     SeqXS2T,SeqXTStar);
%     % 
%     % KTTStar = CovCal(CovFun,thetaTarget,XTStar,XTStar,...
%     %     SeqXTStar,SeqXTStar);
%     % KTTStar = cell2mat(KTTStar);
%     % KTTStar = reshape(KTTStar,NStar,NStar,[]);
%     % S4Minus = sigmaTarget(1)^2*eye(NStar);
%     % KTTStar = KTTStar-S4Minus;
%     % KTTStar = sum(KTTStar,3);
% 
%     KT_KST = KTarget{1}+1e-6*eye(sum(SeqXTarget{1}));
%     L_KT_KST = chol(KT_KST,'lower');
% 
%     Inv_L_KT_KST = L_KT_KST\eye(m);
%     Inv_KT_KST = Inv_L_KT_KST'*Inv_L_KT_KST; %(KTT-KTS/KSS*KST)^-1
% 
%     yKT = Inv_KT_KST*y{end};
%     KTStarT = KTStar{1};
%     Mu = KTStarT'*yKT;
%     KKTTStar = KTStarT'*Inv_KT_KST*KTStarT;
%     Sigma = KTTStar-KKTTStar;
% 
% end
% 
% end
% 
% function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2)
%     K = cellfun(CovFun,theta,X1,X2,...
%             SeqX1,SeqX2,'UniformOutput',false);
% end
% 
% function KK = Prod(Kii,Kij)
% 
% KK = Kij'/Kii;
% 
% end
% 
% function KK = ProdV(G,Kij)
% 
% KK = G*Kij;
% 
% end
% 
% function KK = ProdVV(G,Inv_KT_KST)
% 
% KK = -G'*Inv_KT_KST;
% 
% end

% function yK = yS_KST(yS,Inv_KST)
% yK = Inv_KST'*yS;
% end


% function yK = yT_KST(yT,Inv_KST)
% yK = Inv_KST*yT;
% end