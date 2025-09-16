function [KTarget, KSTInd, KTTInd, Gxx_n_n, Inv_KS, A] =...
    CovCalInd2D(theta,y,Xf,XStar,CovType,NDomain,NStream)

% NDomain = length(y);

CovFun = str2func(CovType);

%% Prepare model parameters
D = size(Xf{1},2);
thetaSource = cell(1,NDomain-1);
thetaST = cell(1,NDomain-1);
thetaTarget = cell(1,NDomain);
switch CovType
    case 'CovCP2D'
        temp = [NStream;repmat(NStream(end),NDomain-1,1)];
        lRho = sum(NStream(1:NDomain-1))+NDomain*NStream(end);
        Rho = theta(1:lRho);
        Rho = Rho(:);
        Rho = mat2cell(Rho,temp);
        lScale = lRho*D;
        Scale = theta(lRho+(1:lScale));
        Scale = Scale(:);
        Scale = mat2cell(Scale,D*temp);
        sigma = theta(lRho+lScale+(1:NDomain));
        for i = 1:NDomain
            if i<NDomain
                thetaSource{i} = [Rho{i};Scale{i};sigma(i)];
                thetaTarget{i} = [Rho{i+NDomain-1};Scale{i+NDomain-1};sigma(end)];
                % thetaS2T{i} = [[Rho{i};Scale{i};sigma(i)] ...
                %                     [Rho{i+NDomain-1};Scale{i+NDomain-1};sigma(i)]];
                thetaST{i} = [Rho{i};Rho{i+NDomain-1};...
                                Scale{i};Scale{i+NDomain-1};sigma(i)];
            else
                thetaTarget{i} = [Rho{i+NDomain-1};Scale{i+NDomain-1};sigma(i)];
            end
        end
        thetaST = [thetaST thetaTarget(end)];
        sigmaTarget = repmat(theta(end),NStream(end),NDomain);

    case 'CovSGPTransfer2D'
        temp = [NStream;repmat(NStream(end),NDomain-1,1)];
        lRho = sum(NStream(1:NDomain-1))+NDomain*NStream(end);
        Rho = theta(1:lRho);
        Rho = Rho(:);
        Rho = mat2cell(Rho,temp);
        lScale = lRho*D;
        Scale = theta(lRho+(1:lScale));
        Scale = Scale(:);
        Scale = mat2cell(Scale,D*temp);
        sigma = theta(lRho+lScale+(1:NDomain));
        for i = 1:NDomain
            if i<NDomain
                thetaSource{i} = [repmat(Rho{i},2,1);repmat(Scale{i},2,1);sigma(i)];
                thetaTarget{i} = [repmat(Rho{i+NDomain-1},2,1);...
                    repmat(Scale{i+NDomain-1},2,1);sigma(end)];
                thetaST{i} = [Rho{i};Rho{i+NDomain-1};...
                                Scale{i};Scale{i+NDomain-1};sigma(i)];
            else
                thetaTarget{i} = [repmat(Rho{i+NDomain-1},2,1);...
                    repmat(Scale{i+NDomain-1},2,1);sigma(i)];
            end
        end
        thetaST = [thetaST thetaTarget(end)];
        sigmaTarget = repmat(theta(end),NStream(end),NDomain);

    case 'CovLMC2D'
        NLatent = (length(theta)-NDomain)/...
            (sum(NStream)+(NDomain-1)*NStream(end)+NDomain*D);
        temp = [NStream;repmat(NStream(end),NDomain-1,1)];
        lRho = (sum(NStream(1:NDomain-1))+NDomain*NStream(end))*NLatent;
        Rho = theta(1:lRho);
        Rho = Rho(:);
        Rho = mat2cell(Rho,temp*NLatent);
        lScale = NDomain*D*NLatent;
        Scale = theta(lRho+(1:lScale));
        Scale = Scale(:);
        Scale = mat2cell(Scale,D*NLatent*ones(1,NDomain));
        sigma = theta(lRho+lScale+(1:NDomain));
        for i = 1:NDomain
            if i< NDomain
                thetaSource{i} = [repmat(Rho{i},2,1);Scale{i};sigma(i)];
                thetaTarget{i} = [repmat(Rho{i+NDomain-1},2,1);...
                    Scale{i};sigma(end)];
                thetaST{i} = [Rho{i};Rho{i+NDomain-1};...
                    Scale{i};sigma(i)];
            else
                thetaTarget{i} = [repmat(Rho{i+NDomain-1},2,1);...
                    Scale{i};sigma(end)];
            end
        end
        thetaST = [thetaST thetaTarget(end)];
        sigmaTarget = repmat(theta(end),NStream(end),NDomain);

    case 'CovNonTransfer2D'
        lRho = NStream;
        Rho = theta(1:lRho);
        Rho = Rho(:);
        Rho = mat2cell(Rho,NStream);
        lScale = lRho*D;
        Scale = theta(lRho+(1:lScale));
        Scale = Scale(:);
        Scale = mat2cell(Scale,D*NStream);
        sigma = theta(lRho+lScale+(1:NDomain));
        for i = 1:NDomain
            thetaTarget{i} = [repmat(Rho{i},2,1);...
                repmat(Scale{i},2,1);sigma(end)];
        end
        thetaST = [thetaST thetaTarget(end)];
        sigmaTarget = repmat(theta(end),NStream(end),NDomain);

    otherwise
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

        thetaSource = [rhoSource;lambdaSource;sigmaSource];
        thetaSource = mat2cell(thetaSource,5*NStream,ones(1,NDomain-1));
        thetaTarget = [rhoTarget;lambdaTarget;sigmaTarget];
        thetaTarget = mat2cell(thetaTarget,5*NStream,ones(1,NDomain));
end

%%

XSource = Xf(1:NDomain-1);
XTarget = Xf(NDomain);
XTInd = XStar;

SeqXSource = cellfun(@length,XSource)./NStream(1:NDomain-1)';
SeqXSource = repelem(SeqXSource(:),NStream(1:NDomain-1),1);
SeqXSource = mat2cell(SeqXSource,NStream(1:NDomain-1))';

SeqXTarget = repmat(length(XTarget{1})./NStream(NDomain),NStream(NDomain),NDomain);
SeqXTarget = mat2cell(SeqXTarget,NStream(NDomain),ones(1,NDomain));

SeqXTStar = repmat(length(XTInd{1})/NStream(NDomain),NStream(NDomain),NDomain);
SeqXTStar = mat2cell(SeqXTStar,NStream(NDomain),ones(1,NDomain));

XS2TInd = [XSource XTInd];
SeqXS2TInd =[SeqXSource SeqXTStar{1}];

XS2Target = [XSource XTarget];
SeqXS2Target =[SeqXSource SeqXTarget{1}];

XTarget = repmat(XTarget,1,NDomain);
XTInd = repmat(XTInd,1,NDomain);

% thetaST = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
%     'UniformOutput',false);
% thetaST = [thetaST thetaTarget(end)];

% \Sigma_xx^{n,n}
KTarget = CovCal(CovFun,thetaST(end),XTarget(end),XTarget(end),...
    SeqXTarget(end),SeqXTarget(end));
KTarget = KTarget{end};

% \Gamma_sx^{n,<n} \Gamma_ss^{n,n}
KSTInd = CovCal(CovFun,thetaST,XS2TInd,XTInd,...
    SeqXS2TInd,SeqXTStar);

% \Gamma_sx^{n,n}
KTTInd = CovCal(CovFun,thetaST(end),XTarget(end),XTInd(end),...
    SeqXTarget(end),SeqXTStar(end));
KTTInd = KTTInd{end};

% \Gamma_xx^{n,<n}
KST = CovCal(CovFun,thetaST(1:end-1),XS2Target(1:end-1),XTarget(1:end-1),...
    SeqXS2Target(1:end-1),SeqXTarget(1:end-1));

[mT,nT] = size(KTarget);
[mTS,nTS] = size(KSTInd{end});
S4Minus = sigmaTarget(1)^2*eye(nTS);
KSTInd{end} = KSTInd{end}-S4Minus;

if strcmp(CovType,'CovCP2D')
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXTarget{end},SeqXTarget{end});
    KTarget = KTarget.*Temp;
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXS2TInd{end},SeqXTStar{end});
    KSTInd{end} = KSTInd{end}.*Temp;
    Temp = eye(NStream(end));
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
        SeqXSource,SeqXSource);
    % inv(\Sigma_xx^{<n,<n})
    Inv_KS = cellfun(@inv,KSource,'UniformOutput',false);

    % \Sigma_xx^{n,n}
    XTT = XTarget(1:end-1);
    thetaTT = thetaTarget(1:NDomain-1);
    SeqXTT = SeqXTarget(1:end-1);
    KTargetAdd = CovCal(CovFun,thetaTT,XTT,XTT,...
        SeqXTT,SeqXTT);
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
        SeqXTTStar,SeqXTTStar);
    KTStarAdd = cell2mat(KTStarAdd);
    KTStarAdd = reshape(KTStarAdd,mTS,nTS,[]);
    S4Minus = sigmaTarget(1)^2*eye(nTS);
    KTStarAdd = KTStarAdd-S4Minus;
    KTStarAdd = sum(KTStarAdd,3);
    KSTInd{end} = KSTInd{end}+KTStarAdd;

    % \Gamma_sx^{n,n}
    KTTStarAdd = CovCal(CovFun,thetaTT,XTarget(1:end-1),XTInd(1:end-1),...
        SeqXTarget(1:end-1),SeqXTStar(1:end-1));
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

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2)
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,'UniformOutput',false);
end

    
