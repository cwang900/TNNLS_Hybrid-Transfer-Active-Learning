function X = OPTIMSE(XSample0,XStar,Sigma,XInd,theta,NDomain,CovType,lb,ub,NStream)

[NS,ND] = size(XSample0);
opt.algorithm = NLOPT_LD_MMA;%NLOPT_LN_AUGLAG_EQ;
opt.maxeval = 10000;
opt.lower_bounds = lb*ones(1,ND);
opt.upper_bounds = ub*ones(1,ND);
opt.xtol_rel = 1e-4;

opt.min_objective = @(XSample) IMSE(XSample,XStar,Sigma,XInd,theta,NDomain,CovType,NStream);
J = zeros(1,NS);
XX = zeros(2,NS);
for i = 1:NS
    XX(:,i) = nlopt_optimize(opt, XSample0(i,:));
    J(i) = IMSE(XX(:,i),XStar,Sigma,XInd,theta,NDomain,CovType,NStream);
end
[~,II] = min(J);
X = XX(:,II);

end

function [J, dJ] = IMSE(XSample,XStar,Sigma,XInd,theta,NDomain,CovType,NStream)

XSample0 = XSample;
theta0 = theta;
XInd0 = XInd;
XStar0 = XStar;
XSample = repmat(XSample(:),NStream,1);

yt = zeros(size(XSample));
alpha = zeros(size(XInd));
% [~, SigmaT] = PostOnline(theta,yt,XSample,Mu,Sigma,XInd,NDomain,CovType);
[~, SigmaT] = PostInd(theta,{yt},{XSample},{XStar},alpha,Sigma,{XInd},NDomain,CovType,NStream);
% J = OptJAL(theta,XInd,XInd,Sigma,NDomain,CovType);

%% Numerical integration

% if strcmp(CovType,'CovOneShare')
%     S = diag(SigmaT);
%     S1 = S(1:end/2);
%     S2 = S(end/2+1:end);
%     S12 = SigmaT(1:end/2,end/2+1:end);
%     S12 = diag(S12);
%     Y = S1.*S2-S12.^2;
%     X = XInd(1:end/2);
%     J = log(trapz(X,Y));
%     SigmaT = SigmaT+1e-6*eye(size(SigmaT,1));
%     SL = chol(SigmaT,'lower');
%     J = 2*sum(log(diag(SL)));
%     J = det(SigmaT);
% else
%     X1 = XInd(1:end/2);
%     X2 = XInd(end/2+1:end);
%     S = diag(SigmaT);
%     Y1 = S(1:end/2);
%     Y2 = S(end/2+1:end);
%     J = trapz(X1,Y1)+trapz(X2,Y2);
% end
% X1 = XStar(1:end/2);
% X2 = XStar(end/2+1:end);
% S = diag(SigmaT);
% Y1 = S(1:end/2);
% Y2 = S(end/2+1:end);
% J = trapz(X1,Y1)+trapz(X2,Y2);
XStar = reshape(XStar,[],NStream);
S = diag(SigmaT);
S = reshape(S,[],NStream);
J= 0;
for i = 1:NStream
    XX = XStar(:,i);
    YY = S(:,i);
    J = J+trapz(XX,YY);
end
J = log(J);

%%

% J = sum(diag(SigmaT));

if nargout == 2

    fun = @(XSample)IMSE(XSample,XStar0,Sigma,XInd0,theta0,NDomain,CovType,NStream);
    dJ = Grad(fun, XSample0);

end

end


%% Debug part

% function [J, dJ] = IMSE(XSample,XStar,X,Xu,theta,CovType)
% 
% XSample0 = XSample;
% theta0 = theta;
% NDomain = length(X);
% % XStar = [XStar;XStar];
% NStar = length(XStar);
% XSample = XSample(:);
% X{end} = reshape(X{end},[],2);
% X{end} = [X{end}; [XSample XSample]];
% X{end} = X{end}(:);
% 
% %% Prepare model parameters
% 
% switch CovType
%     case 'CovOneShare'
%         rhoSource = reshape(theta(1:2*(NDomain-1)),2,[]);
%         rhoSource = repelem(rhoSource,2,1);
%         rhoSource([2 4],:) = 0;
%         theta(1:2*(NDomain-1)) = [];
% 
%         rhoTarget = reshape(theta(1:2*(NDomain)),2,[]);
%         rhoTarget = repelem(rhoTarget,2,1);
%         rhoTarget([2 4],:) = 0;
%         theta(1:2*(NDomain)) = [];
% 
%         lambdaSource = reshape(theta(1:2*(NDomain-1)),2,[]);
%         lambdaSource = repelem(lambdaSource,2,1);
%         theta(1:2*(NDomain-1)) = [];
% 
%         lambdaTarget = reshape(theta(1:2*(NDomain)),2,[]);
%         lambdaTarget = repelem(lambdaTarget,2,1);
%         theta(1:2*(NDomain)) = [];
% 
%         sigmaSource = repmat(theta(1:1*(NDomain-1))',2,1);
%         theta(1:1*(NDomain-1)) = [];
%         sigmaTarget = repmat(theta,2,NDomain);
%     case 'CovSGPTransfer'
%         rhoSource = reshape(theta(1:2*(NDomain-1)),2,[]);
%         rhoSource = repelem(rhoSource,2,1);
%         rhoSource([2 4],:) = 0;
%         theta(1:2*(NDomain-1)) = [];
% 
%         rhoTarget = reshape(theta(1:2*(NDomain)),2,[]);
%         rhoTarget = repelem(rhoTarget,2,1);
%         rhoTarget([2 4],:) = 0;
%         theta(1:2*(NDomain)) = [];
% 
%         lambdaSource = reshape(theta(1:2*(NDomain-1)),2,[]);
%         lambdaSource = repelem(lambdaSource,2,1);
%         theta(1:2*(NDomain-1)) = [];
% 
%         lambdaTarget = reshape(theta(1:2*(NDomain)),2,[]);
%         lambdaTarget = repelem(lambdaTarget,2,1);
%         theta(1:2*(NDomain)) = [];
% 
%         sigmaSource = repmat(theta(1:1*(NDomain-1))',2,1);
%         theta(1:1*(NDomain-1)) = [];
%         sigmaTarget = repmat(theta,2,NDomain);
% 
%     case 'CovNonTransfer'
%         theta = theta(:);
%         rho = repelem(theta(1:2),2,1);
%         rho([2 4],:) = 0;
%         lambda = repelem(theta(3:4),2,1);
%         lambda([2 4],:) = 0;
%         sigma = repmat(theta(end),2,1);
%         theta = [rho;lambda;sigma];
% end
% 
% if ~strcmp(CovType,'CovNonTransfer')
% 
% thetaSource = [rhoSource;lambdaSource;sigmaSource];
% thetaSource = mat2cell(thetaSource,10,ones(1,NDomain-1));
% thetaTarget = [rhoTarget;lambdaTarget;sigmaTarget];
% thetaTarget = mat2cell(thetaTarget,10,ones(1,NDomain));
% 
% %% Calculate each covariance block
% 
% fun = str2func(CovType);
% 
% XSource = X(1:NDomain-1);
% SeqXSource = cellfun(@length,XSource)/2;
% SeqXSource = repmat(SeqXSource,2,1);
% SeqXSource = mat2cell(SeqXSource,2,ones(1,NDomain-1));
% XTarget = X{NDomain};
% SeqXTarget = repmat(length(XTarget)/2,2,NDomain);
% SeqXTarget = mat2cell(SeqXTarget,2,ones(1,NDomain));
% KSource = cellfun(fun,thetaSource,XSource,XSource,...
%     SeqXSource,SeqXSource,'UniformOutput',false);
% 
% XS2T = [XSource {XTarget}];
% SeqXS2T =[SeqXSource SeqXTarget{1}];
% XTarget = repmat(XTarget,1,NDomain);
% XTarget = mat2cell(XTarget,size(XTarget,1),ones(1,NDomain));
% thetaST = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
%     'UniformOutput',false);
% thetaST = [thetaST thetaTarget(end)];
% KTarget = cellfun(fun,thetaST,XS2T,XTarget,...
%     SeqXS2T,SeqXTarget,'UniformOutput',false);
% [m,n] = size(KTarget{end});
% 
% XTT = XTarget(1:end-1);
% thetaTT = thetaTarget(1:NDomain-1);
% SeqXTT = SeqXTarget(1:end-1);
% KTargetAdditional = cellfun(fun,thetaTT,XTT,XTT,...
%     SeqXTT,SeqXTT,'UniformOutput',false);
% KTargetAdditional = cell2mat(KTargetAdditional);
% KTargetAdditional = reshape(KTargetAdditional,m,n,[]);
% S4Minus = sigmaTarget(1)^2*eye(m);
% KTargetAdditional = KTargetAdditional-S4Minus;
% KTargetAdditional = sum(KTargetAdditional,3);
% KTarget{end} = KTarget{end}+KTargetAdditional;
% 
% XTStar = repmat({XStar},1,NDomain);
% thetaTStar = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
%     'UniformOutput',false);
% thetaTStar = [thetaTStar, thetaTarget(end)];
% SeqXTStar = repmat(NStar/2,2,NDomain);
% SeqXTStar = mat2cell(SeqXTStar,2,ones(1,NDomain));
% KTStar = cellfun(fun,thetaTStar,X,XTStar,...
%     SeqXS2T,SeqXTStar,'UniformOutput',false);
% KTStarAdditional = cellfun(fun,thetaTarget(1:NDomain-1),...
%     XTT,XTStar(1:NDomain-1),...
%     SeqXTT,SeqXTStar(1:NDomain-1),'UniformOutput',false);
% KTStarAdditional = cell2mat(KTStarAdditional);
% KTStarAdditional = reshape(KTStarAdditional,sum(SeqXTT{1}),sum(SeqXTStar{1}),[]);
% KTStarAdditional = sum(KTStarAdditional,3);
% KTStar{end} = KTStar{end}+KTStarAdditional;
% 
% KTTStar = cellfun(fun,thetaTarget,XTStar,XTStar,...
%     SeqXTStar,SeqXTStar,'UniformOutput',false);
% KTTStar = cell2mat(KTTStar);
% KTTStar = reshape(KTTStar,NStar,NStar,[]);
% S4Minus = sigmaTarget(1)^2*eye(NStar);
% KTTStar = KTTStar-S4Minus;
% KTTStar = sum(KTTStar,3);
% 
% KTStarT = KTStar{end};
% 
% %% 
% 
% G = cellfun(@Prod,KSource,KTarget(1:NDomain-1),'UniformOutput',false);
% Gtemp = cell2mat(G);
% % KSKTKS = Gtemp*cell2mat(KTarget(1:NDomain-1));
% KST = cellfun(@ProdV,G,KTarget(1:NDomain-1),'UniformOutput',false);
% KST = cell2mat(KST);
% KST = reshape(KST,m,n,[]);
% KST = sum(KST,3);
% % temp = KTargetAdditional-KST;
% % if temp(1) < 0
% %     aa;
% % end
% KT_KST = KTarget{end}-KST;
% KT_KST = KT_KST+1e-6*eye(m);
% L_KT_KST = chol(KT_KST,'lower');
% 
% Inv_L_KT_KST = L_KT_KST\eye(m);
% Inv_KT_KST = Inv_L_KT_KST'*Inv_L_KT_KST; %(KTT-KTS/KSS*KST)^-1
% 
% Inv_KS = cellfun(@inv,KSource,'UniformOutput',false);
% KST_InvOT = cellfun(@(G)ProdVV(G,Inv_KT_KST),G,...
%     'UniformOutput',false);
% % yKST_InvOT = cellfun(@(G)ProdVV(G,Inv_KT_KST),G,...
% %     'UniformOutput',false);
% KST_InvOT_mat = cell2mat(KST_InvOT');
% Inv_KS_Sub = -KST_InvOT_mat*Gtemp;
% Inv_KS_Full = blkdiag(Inv_KS{:})+Inv_KS_Sub;
% 
% %%
% 
% KTStarS = cell2mat(KTStar(1:NDomain-1)');
% K11 = KTStarS'*Inv_KS_Full*KTStarS;
% K12 = KTStarT'*KST_InvOT_mat'*KTStarS;
% K22 = KTStarT'*Inv_KT_KST*KTStarT;
% KAll = KTTStar-(K11+K12+K12'+K22);
% 
% else
% 
%     SeqX = length(X{end})/2;
%     SeqX = [SeqX;SeqX];
% K = CovOneShare(theta,X{end},X{end},SeqX,SeqX);
% 
% SeqXS = length(XStar)/2;
% SeqXS = [SeqXS;SeqXS];
% KS = CovOneShare(theta,X{end},XStar,SeqX,SeqXS);
% KSS = CovOneShare(theta,XStar,XStar,SeqXS,SeqXS);
% 
% KAll = KSS-KS'/K*KS;
% 
% end
% 
% J = sum(diag(KAll));
% 
% if nargout == 2
% 
%     fun = @(XSample)IMSE(XSample,XStar,X,Xu,theta0,CovType);
%     dJ = Grad(fun, XSample0);
% 
% end
% 
% end


function KK = Prod(Kii,Kij)

KK = Kij'/Kii;

end

function KK = ProdV(G,Kij)

KK = G*Kij;

end

function KK = ProdVV(G,Inv_KT_KST)

KK = -G'*Inv_KT_KST;

end

%% Prepare input for covariance function

% rho = reshape(theta(1:2*(NDomain)),2,[]);
% rho = repelem(rho,2,1);
% rho([2 4],:) = 0;
% 
% theta(1:2*(NDomain)) = [];
% 
% lambda = reshape(theta(1:2*(NDomain)),2,[]);
% lambda = repelem(lambda,2,1);
% theta(1:2*(NDomain)) = [];
% 
% sigma = repmat(theta,2,NDomain);
% 
% theta = [rho;lambda;sigma];
% theta = mat2cell(theta,10,ones(1,NDomain));
% 
% SeqXt = repmat(length(XSample)/2,2,NDomain);
% SeqXt = mat2cell(SeqXt,2,ones(1,NDomain));
% Nt = length(XSample);
% XSample = repmat({XSample},1,NDomain);
% 
% SeqXInd = repmat(length(XInd)/2,2,NDomain);
% SeqXInd = mat2cell(SeqXInd,2,ones(1,NDomain));
% NInd = length(XInd);
% XInd = repmat({XInd},1,NDomain);
% 
% %% Calculate the covariance functions
% 
% if strcmp(CovType,'CovSGPTransfer') 
%         fun = str2func(CovType);
%     else
%         fun = str2func('CovOneShare');
% end
% 
% KTemp = cellfun(fun,theta,XSample,XSample,...
%     SeqXt,SeqXt,'UniformOutput',false);
% KTemp = cell2mat(KTemp);
% KTemp = reshape(KTemp,Nt,Nt,[]);
% % S4Minus = sigma(1)^2*eye(Nt);
% % KTemp = KTemp-S4Minus;
% Ktt = sum(KTemp,3);
% 
% KTemp = cellfun(fun,theta,XInd,XInd,...
%     SeqXInd,SeqXInd,'UniformOutput',false);
% KTemp = cell2mat(KTemp);
% KTemp = reshape(KTemp,NInd,NInd,[]);
% S4Minus = sigma(1)^2*eye(NInd);
% KTemp = KTemp-S4Minus;
% KXX = sum(KTemp,3)+1e-6*eye(NInd);
% 
% KTemp = cellfun(fun,theta,XSample,XInd,...
%     SeqXt,SeqXInd,'UniformOutput',false);
% KTemp = cell2mat(KTemp);
% KTemp = reshape(KTemp,Nt,NInd,[]);
% KXtXInd = sum(KTemp,3);
% 
% Gt = KXtXInd/KXX;
% Pt = Ktt+Gt*(Sigma-KXX)*Gt';
% PtL = chol(Pt,"upper");
% SGP = Sigma*Gt'/PtL;
% 
% SigmaT = Sigma-SGP*SGP';
% J = sum(diag(SigmaT));