function [Gxs_nn, Qxx_nn, KTStar] =...
    CovCalPred2D(theta,XInd,XStar,NDomain,CovType,NStream)

CovFun = str2func(CovType);

%% Prepare model parameters
D = size(XInd{1},2);
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

XTInd = XInd;
XTStar = XStar;

SeqXInd = repmat(length(XTInd{1})/NStream(NDomain),NStream(NDomain),NDomain);
SeqXInd = mat2cell(SeqXInd,NStream(NDomain),ones(1,NDomain));

SeqXTStar = repmat(length(XTStar{1})/NStream(NDomain),NStream(NDomain),NDomain);
SeqXTStar = mat2cell(SeqXTStar,NStream(NDomain),ones(1,NDomain));

XTInd = repmat(XTInd,1,NDomain);
XTStar = repmat(XTStar,1,NDomain);

% \Gamma_ss_nn
KTInd = CovCal(CovFun,thetaTarget(end),XTInd(end),XTInd(end),...
    SeqXInd(end),SeqXInd(end));
KTInd = KTInd{end};

% \Gamma_sx_nn
KTIndTStar = CovCal(CovFun,thetaTarget(end),XTInd(end),XTStar(end),...
    SeqXInd(end),SeqXTStar(end));
KTIndTStar = KTIndTStar{end};

% \Gamma_xx_nn
KTStar = CovCal(CovFun,thetaTarget(end),XTStar(end),XTStar(end),...
    SeqXTStar(end),SeqXTStar(end));
KTStar = KTStar{end};

[mT,nT] = size(KTInd);
[mTS,nTS] = size(KTStar);
S4Minus = sigmaTarget(1)^2*eye(mT);
KTInd = KTInd-S4Minus;
if strcmp(CovType,'CovCP2D')
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXInd{end},SeqXInd{end});
    KTInd = KTInd.*Temp;
    Temp = eye(NStream(end));
    Temp = repelem(Temp,SeqXTStar{end},SeqXTStar{end});
    KTStar = KTStar.*Temp;
    Temp = eye(NStream(end));
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
        SeqXTT,SeqXTT);
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
        SeqXTTStar,SeqXTTStar);
    KTStarAdd = cell2mat(KTStarAdd);
    KTStarAdd = reshape(KTStarAdd,mTS,nTS,[]);
    S4Minus = sigmaTarget(1)^2*eye(nTS);
    KTStarAdd = KTStarAdd-S4Minus;
    KTStarAdd = sum(KTStarAdd,3);
    KTStar = KTStar+KTStarAdd-S4Minus;

    % \Gamma_sx_nn
    KTIndTStarAdd = CovCal(CovFun,thetaTT,XTInd(1:end-1),XTStar(1:end-1),...
        SeqXInd(1:end-1),SeqXTStar(1:end-1));
    KTIndTStarAdd = cell2mat(KTIndTStarAdd);
    KTIndTStarAdd = reshape(KTIndTStarAdd,mT,mTS,[]);
    KTIndTStarAdd = sum(KTIndTStarAdd,3);
    KTIndTStar = KTIndTStar+KTIndTStarAdd;
    
end
KTInd = KTInd+1e-6*eye(mT);
Gxs_nn = KTIndTStar'/KTInd;
Qxx_nn = Gxs_nn*KTIndTStar;

end

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2)
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,'UniformOutput',false);
end


    
