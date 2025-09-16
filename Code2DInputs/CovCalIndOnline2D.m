function [Gxs_n_n, A] =...
    CovCalIndOnline2D(theta,y,Xf,XStar,CovType,NDomain,NStream)

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

XTarget = Xf(end);
XTInd = XStar;

SeqXTarget = repmat(length(XTarget{1})/NStream(NDomain),NStream(NDomain),NDomain);
SeqXTarget = mat2cell(SeqXTarget,NStream(NDomain),ones(1,NDomain));

SeqXTInd = repmat(length(XTInd{1})/NStream(NDomain),NStream(NDomain),NDomain);
SeqXTInd = mat2cell(SeqXTInd,NStream(NDomain),ones(1,NDomain));

XTarget = repmat(XTarget,1,NDomain);
XTInd = repmat(XTInd,1,NDomain);

% \Sigma_xx^{n,n}
KTarget = CovCal(CovFun,thetaTarget(end),XTarget(end),XTarget(end),...
    SeqXTarget(end),SeqXTarget(end));
KTarget = KTarget{end};

% \Gamma_xx^{n,n}
KTTInd = CovCal(CovFun,thetaTarget(end),XTarget(end),XTInd(end),...
    SeqXTarget(end),SeqXTInd(end));
KTTInd = KTTInd{end};

% \Gamma_sx^{n,n}
KTInd = CovCal(CovFun,thetaTarget(end),XTInd(end),XTInd(end),...
    SeqXTInd(end),SeqXTInd(end));
KTInd = KTInd{end};

[mT,nT] = size(KTarget);
[mTInd,nTInd] = size(KTInd);
S4Minus = sigmaTarget(1)^2*eye(mTInd);
KTInd = KTInd-S4Minus;

if strcmp(CovType,'CovCP2D')
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXTarget{end},SeqXTarget{end});
    KTarget = KTarget.*Temp;
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXTarget{end},SeqXTInd{end});
    KTTInd = KTTInd.*Temp;
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXTInd{end},SeqXTInd{end});
    KTInd = KTInd.*Temp;
end

%%

% For transfer learning based methods
if ~isempty(thetaSource)

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

    % \Gamma_sx^{n,n}
    KTTStarAdd = CovCal(CovFun,thetaTT,XTarget(1:end-1),XTInd(1:end-1),...
        SeqXTarget(1:end-1),SeqXTInd(1:end-1));
    KTTStarAdd = cell2mat(KTTStarAdd);
    KTTStarAdd = reshape(KTTStarAdd,mT,mTInd,[]);
    KTTStarAdd = sum(KTTStarAdd,3);
    KTTInd = KTTInd+KTTStarAdd;

    % \Gamma_ss^{n,n}
    KTIndAdd = CovCal(CovFun,thetaTT,XTInd(1:end-1),XTInd(1:end-1),...
        SeqXTInd(1:end-1),SeqXTInd(1:end-1));
    KTIndAdd = cell2mat(KTIndAdd);
    KTIndAdd = reshape(KTIndAdd,mTInd,nTInd,[]);
    S4Minus = sigmaTarget(1)^2*eye(mTInd);
    KTIndAdd = KTIndAdd-S4Minus;
    KTIndAdd = sum(KTIndAdd,3);
    KTInd = KTInd+KTIndAdd;
    
end

% \G_xs^{n,n}
Gxs_n_n = KTTInd/(KTInd+1e-6*eye(mTInd));

% \Q_xx^{n,n}
Qxx_nn = Gxs_n_n*KTTInd';

% \Sigma_xx^{n,n}-Q_xx^{n,n}
A = KTarget-Qxx_nn;

end

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2)
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,'UniformOutput',false);
end

    
