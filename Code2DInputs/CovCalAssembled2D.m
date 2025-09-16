function [Myn, Inv_Cxx_nn, Inv_KS, Det] =...
    CovCalAssembled2D(theta,y,Xf,CovType,Aim,NStream)

NDomain = length(y);

CovFun = str2func(CovType);

%% Prepare model parameters
D = size(Xf{1},2);
thetaSource = cell(1,NDomain-1);
thetaS2T = cell(1,NDomain-1);
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
                thetaS2T{i} = [Rho{i};Rho{i+NDomain-1};...
                                Scale{i};Scale{i+NDomain-1};sigma(i)];
            else
                thetaTarget{i} = [Rho{i+NDomain-1};Scale{i+NDomain-1};sigma(i)];
            end
        end
        thetaS2T = [thetaS2T thetaTarget(end)];
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
                thetaS2T{i} = [Rho{i};Rho{i+NDomain-1};...
                                Scale{i};Scale{i+NDomain-1};sigma(i)];
            else
                thetaTarget{i} = [repmat(Rho{i+NDomain-1},2,1);...
                    repmat(Scale{i+NDomain-1},2,1);sigma(i)];
            end
        end
        thetaS2T = [thetaS2T thetaTarget(end)];
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
                thetaS2T{i} = [Rho{i};Rho{i+NDomain-1};...
                    Scale{i};sigma(i)];
            else
                thetaTarget{i} = [repmat(Rho{i+NDomain-1},2,1);...
                    Scale{i};sigma(end)];
            end
        end
        thetaS2T = [thetaS2T thetaTarget(end)];
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
        thetaS2T = [thetaS2T thetaTarget(end)];
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

SeqXSource = cellfun(@length,XSource)./NStream(1:NDomain-1)';
SeqXSource = repelem(SeqXSource(:),NStream(1:NDomain-1),1);
SeqXSource = mat2cell(SeqXSource,NStream(1:NDomain-1))';

SeqXTarget = repmat(length(XTarget{1})./NStream(NDomain),NStream(NDomain),NDomain);
SeqXTarget = mat2cell(SeqXTarget,NStream(NDomain),ones(1,NDomain));

XS2T = [XSource XTarget];
SeqXS2T =[SeqXSource SeqXTarget{1}];
XTarget = repmat(XTarget,1,NDomain);

% thetaS2T = cellfun(@horzcat,thetaSource,thetaTarget(1:NDomain-1),...
%     'UniformOutput',false);
% thetaS2T = [thetaS2T thetaTarget(end)];

% \Sigma_xx^{n,n}
KTarget = CovCal(CovFun,thetaS2T(end),XTarget(end),XTarget(end),...
    SeqXTarget(end),SeqXTarget(end));
KTarget = KTarget{end};

% \Gamma_xx^{<n,n}
KS2T = CovCal(CovFun,thetaS2T,XS2T,XTarget,...
    SeqXS2T,SeqXTarget);

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
    L_Cxx_nn = jitChol(Cxx_nn)';%chol(Cxx_nn,'lower');
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

function K = CovCal(CovFun,theta,X1,X2,SeqX1,SeqX2)
    NStream = [length(SeqX1); length(SeqX2)];
    NDomain = length(theta);
    NStream = repmat({NStream},1,NDomain);
    K = cellfun(CovFun,theta,X1,X2,...
            SeqX1,SeqX2,'UniformOutput',false);
end

    
